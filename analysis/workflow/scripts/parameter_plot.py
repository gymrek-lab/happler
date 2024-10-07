#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
import matplotlib
import numpy as np
matplotlib.use('Agg')
import numpy.typing as npt
from haptools import logging
import matplotlib.pyplot as plt
from haptools.ld import pearson_corr_ld
from scipy.optimize import linear_sum_assignment
from haptools.data import Genotypes, GenotypesVCF, GenotypesPLINK, Haplotypes

from snakemake_io import glob_wildcards


DTYPES = {
    "beta": np.float64,
    "ld": np.float64,
    "rep": np.uint8,
    "gt": "U5",
    "alpha": np.float64,
    "samp": np.uint32,
    "num_haps": np.uint8,
}

LOG_SCALE = {"alpha",}


def match_haps(gts: Genotypes, observed: Haplotypes, causal: Haplotypes, drop_extras: bool = True) -> tuple:
    """
    Match the observed and causal haplotypes to maximize best pairwise LD

    This function uses scipy.optimize.linear_sum_assignment to find the optimal
    match between causal and observed haplotypes to maximize the sum of the pairwise LD

    Parameters
    ----------
    gts: Genotypes
        The set of genotypes for each of the variants in the haplotypes
    observed: Path
        A path to a .hap file containing haplotypes output by happler
    causal: Path
        A path to a .hap file containing a simulated causal haplotype
    drop_extras: bool, optional
        Drop any observed haplotypes that don't have a match in the causal haps.
        Otherwise, include a match for every observed haplotype!

    Returns
    -------
    npt.NDArray
        An array of the LD values for each pair of observed and causal haplotypes
    npt.NDArray
        An array of the row indices indicating the optimal match of causal haplotypes
        to observed haplotypes
    npt.NDArray
        An array of the column indices indicating the optimal match of causal haplotypes
        to observed haplotypes
    """
    obs = observed.transform(gts).data.sum(axis=2)
    exp = causal.transform(gts).data.sum(axis=2)
    # Compute the pairwise LD between each observed and causal haplotype
    # The rows and columns of ld_mat both correspond to the observed and causal haps in
    # the order they were given. For example, if there is 1 observed hap and 3 causal
    # haps, then there should be four rows and four columns in ld_mat
    # Note that ld_mat is symmetric, so the upper triangle is always redundant
    # TODO: use the latest haptools.ld.pearson_corr_ld to do this, instead
    ld_mat = np.corrcoef(obs, exp, rowvar=False)
    # Retrieve only the correlation of observed haps (as rows) vs causal haps (as cols)
    ld_mat = np.abs(ld_mat[:obs.shape[1], obs.shape[1]:])
    # Return the ld_mat idxs of the optimal match of each observed hap to a causal hap
    row_idx, col_idx = linear_sum_assignment(ld_mat, maximize=True)
    if not drop_extras:
        missing_rows = np.setdiff1d(range(ld_mat.shape[0]), row_idx)
        # If there are more observed haps than causal haps, then we just assign the
        # remaining haps to the best causal hap that matches for each of them
        row_idx = np.concatenate((row_idx, missing_rows))
        col_idx = np.concatenate((col_idx, ld_mat[missing_rows].argmax(axis=1)))
    return ld_mat, row_idx, col_idx


def get_best_ld(
        gts: GenotypesVCF,
        observed_hap: Path,
        causal_hap: Path,
        region: str = None,
        observed_id: str = None,
        causal_id: str = None,
        show_extras: bool = False,
        log: Logger = None
    ):
    """
    Compute the best LD between the observed haplotypes and the causal haplotype

    Parameters
    ----------
    gts: GenotypesVCF
        A GenotypesVCF (or GenotypesPLINK) object from which to load genotypes
    observed_hap: Path
        A path to a .hap file containing haplotypes output by happler
    causal_hap: Path
        A path to a .hap file containing a simulated causal haplotype
    region: str, optional
        Only load genotypes from this region (or the entire file, otherwise)
    observed_id: str, optional
        The ID to load from the observed_hap file (or just all of the IDs, otherwise)
    causal_id: str, optional
        The ID to load from the causal_hap file (or just the first hap, otherwise)
    show_extras: bool, optional
        Whether to show observed haps that don't match with a causal hap
    log: Logger, optional
        A logging module to pass to haptools
    """
    # load the causal haplotype given by 'causal_id' or just the first hap
    causal_hap = Haplotypes(causal_hap, log=log)
    causal_hap.read(haplotypes=(set((causal_id,)) if causal_id is not None else None))

    # load the observed haplotype given by 'observed_id' or just all of the haps
    observed_hap = Haplotypes(observed_hap, log=log)
    observed_hap.read(haplotypes=(set((observed_id,)) if observed_id is not None else None))

    # load the genotypes
    variants = {
        v.id for hap in causal_hap.data
        for v in causal_hap.data[hap].variants
    }
    variants.update({
        v.id for hap in observed_hap.data
        for v in observed_hap.data[hap].variants
    })
    gts._var_idx = None
    gts.read(variants=variants, region=region)
    gts.check_phase()
    gts.check_missing()
    gts.check_biallelic()

    # compute LD between every observed haplotype and the causal haplotype
    observed_ld, best_row_idx, best_col_idx = match_haps(
        gts, observed_hap, causal_hap, drop_extras=(not show_extras),
    )
    # return the strongest possible LD for each observed hap with a causal hap
    return observed_ld[best_row_idx, best_col_idx]


def get_finemap_metrics(metrics_path: Path, log: Logger = None):
    """
    Parse data from a TSV containing metrics from fine-mapping

    Parameters
    ----------
    metrics: Path
        The path to a metrics.tsv file
    log: Logger, optional
        A logging module to pass to haptools
    """
    # The metrics are:
    # 1) What is the PIP of the observed hap?
    # 2) Does the observed hap get the highest PIP?
    # 3) What is the best PIP among the variants?
    # 4) Is the observed hap in a credible set?
    # 5) What is the purity of the credible set?
    return np.loadtxt(
        fname=metrics_path,
        delimiter=" ",
        dtype=[
            ("pip", np.float64),
            ("has_highest_pip", np.bool_),
            ("best_variant_pip", np.float64),
            ("in_credible_set", np.bool_),
            ("cs_purity", np.float64),
        ],
    )


def count_shared(observed, causal, log, observed_id: str = None, causal_id: str = None):
    observed = Haplotypes(observed, log=log)
    observed.read(haplotypes=(set((observed_id,)) if observed_id is not None else None))
    obs_vars = set(
        var.id
        for hap in observed.data
        for var in observed.data[hap].variants
    )

    causal = Haplotypes(causal, log=log)
    causal.read(haplotypes=(set((causal_id,)) if causal_id is not None else None))
    causal_vars = set(
        var.id
        for hap in causal.data
        for var in causal.data[hap].variants
    )

    # return the number of shared SNPs and the number of expected SNPs
    return len(obs_vars & causal_vars), len(causal_vars)


def plot_params(
        params: npt.NDArray,
        vals: npt.NDArray,
        val_title: str,
    ):
    """
    Plot vals against parameter values

    Parameters
    ----------
    params: npt.NDArray
        A numpy array of mixed dtype. Each column contains the values of a parameter
    vals: npt.NDArray
        A numpy array containing the values of the plot
    val_title: str
        The name of the value that we are plotting
    """
    figsize = matplotlib.rcParams["figure.figsize"]
    if len(params) > 75:
        figsize[0] = len(params) / 15
    if len(params.dtype) > 5:
        figsize[1] = len(params.dtype) * 1.25
    fig, axs = plt.subplots(
        nrows=len(params.dtype)+1, ncols=1,
        sharex=True, figsize=figsize,
    )
    fig.subplots_adjust(hspace=0)
    # create a plot for the vals, first
    x_vals = [j for j, arr in enumerate(vals) for i in arr]
    axs[0].plot(x_vals, np.concatenate(vals), "go", markersize=2)
    axs[0].set_ylabel(val_title, rotation="horizontal", ha="right")
    axs[0].set_xticks(range(0, len(vals)))
    axs[0].set_xticklabels([])
    axs[0].grid(axis="x")
    axs[0].set_ylim(None, 1)
    # now, plot each of the parameter values on the other axes
    for idx, param in enumerate(params.dtype.names):
        val_title = param
        if len(np.unique(params[param])) == 1:
            axs[idx+1].remove()
            continue
        if param in LOG_SCALE:
            params[param] = -np.log10(params[param])
            val_title = "-log " + val_title
        axs[idx+1].plot(params[param], "-")
        axs[idx+1].set_ylabel(val_title, rotation="horizontal", ha="right")
        axs[idx+1].grid(axis="x")
    return fig


def plot_params_simple(
        params: npt.NDArray,
        vals: npt.NDArray,
        val_title: str,
    ):
    """
    Plot vals against only a single column of parameter values

    Parameters
    ----------
    params: npt.NDArray
        A numpy array containing the values of a parameter
    vals: npt.NDArray
        A numpy array containing the values of the plot
    val_title: str
        The name of the value that we are plotting
    """
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True)
    val_xlabel = params.dtype.names[0]
    if val_xlabel in LOG_SCALE:
        params = -np.log10(params)
        val_xlabel = "-log " + val_xlabel
    params = [j for j, arr in zip(params, vals) for i in arr]
    # plot params against vals
    ax.plot(params, np.concatenate(vals), "go", markersize=2)
    ax.set_ylabel(val_title, color="g")
    ax.set_xlabel(val_xlabel)
    return fig


@click.command()
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("observed_hap", type=click.Path(path_type=Path))
@click.argument("causal_hap", type=click.Path(path_type=Path))
@click.option(
    "-m",
    "--metrics",
    type=Path,
    default=None,
    show_default=True,
    help="Plot fine-mapping metrics, as well",
)
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all genotypes",
    help="""
    The region from which to extract genotypes; ex: 'chr1:1234-34566' or 'chr7'\n
    For this to work, the VCF must be indexed and the seqname must match!""",
)
@click.option(
    "-i",
    "--observed-id",
    type=str,
    show_default="the best haplotype in the file",
    help="A haplotype ID from the .hap file output by happler",
)
@click.option(
    "-c",
    "--causal-id",
    type=str,
    show_default="the first haplotype in the file",
    help="A haplotype ID from the causal .hap file",
)
@click.option(
    "--show-extras",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to show observed haps that don't match with a causal hap",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A PNG file containing the desired heatmap",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="ERROR",
    show_default=True,
    help="The level of verbosity desired",
)
def main(
    genotypes: Path,
    observed_hap: Path,
    causal_hap: Path,
    metrics: Path,
    region: str = None,
    observed_id: str = None,
    causal_id: str = None,
    show_extras: bool = False,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "ERROR",
):
    """
    Create a plot summarizing the LD between the observed and causal haplotypes across
    a range of parameters

    Parameters are inferred from brace expressions injected into the paths to the files
    For example, a path like "out/19_45401409-46401409/happler/hap/{beta}/happler.hap"
    will infer the beta values from the portion of the path labeled '{beta}'
    """
    log = logging.getLogger("parameter_plot", level=verbosity)

    # initialize the genotypes class
    gts = GenotypesVCF
    if genotypes.suffix == ".pgen":
        gts = GenotypesPLINK
    gts = gts(genotypes, log=log)

    # extract parameters and parameter values by globbing wildcards
    # params will be a dictionary mapping parameter names to lists of values
    params = dict(glob_wildcards(observed_hap)._asdict())
    # convert the dictionary to a numpy mixed dtype array
    dtypes = {k: DTYPES[k] for k in params.keys()}
    params = np.array(list(zip(*params.values())), dtype=list(dtypes.items()))
    params.sort()

    get_hap_fname = lambda hap_path, param_set: Path(str(hap_path).format(**dict(zip(dtypes.keys(), param_set))))

    # compute LD between the causal hap and the best observed hap across param vals
    ld_vals = [
        get_best_ld(
            gts,
            get_hap_fname(observed_hap, params[idx]),
            get_hap_fname(causal_hap, params[idx]),
            region=region,
            observed_id=observed_id,
            causal_id=causal_id,
            show_extras=show_extras,
            log=log
        )
        for idx in range(len(params))
    ]

    # extract fine-mapping metrics for the observed hap
    if metrics is not None:
        metrics = np.array([
            get_finemap_metrics(
                get_hap_fname(metrics, params[idx]),
                log=log
            )
            for idx in range(len(params))
        ])

    if len(dtypes) > 1 or metrics is not None:
        merged = params
        if metrics is not None:
            merged = np.lib.recfunctions.merge_arrays([params, metrics], flatten=True)
        fig = plot_params(merged, ld_vals, "causal LD")
    elif len(dtypes) == 1:
        fig = plot_params_simple(params, ld_vals, "causal LD")
    else:
        raise ValueError("No parameter values found")

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight", pad_inches=0.05)

if __name__ == "__main__":
    main()
