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
from haptools.data import GenotypesVCF, GenotypesPLINK, Haplotypes

from snakemake_io import glob_wildcards


DTYPES = {
    "beta": np.float64,
    "ld": np.float64,
}


def get_num_haps(hap: Path, log: Logger = None):
    hap = Haplotypes(hap, log=log)
    hap.read()
    return len(hap.data)

def get_best_ld(gts: GenotypesVCF, observed_hap: Path, causal_hap: Path, region: str = None, observed_id: str = None, causal_id: str = None, log: Logger = None):
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
    log: Logger, optional
        A logging module to pass to haptools
    """
    # load the causal haplotype given by 'causal_id' or just the first hap
    causal_hap = Haplotypes(causal_hap, log=log)
    causal_hap.read(haplotypes=(set((causal_id,)) if causal_id is not None else None))
    if causal_id is None:
        causal_id = list(causal_hap.data.keys())[0]
        causal_hap.subset(haplotypes=(causal_id,), inplace=True)
    causal_hap = causal_hap.data[causal_id]

    # load the observed haplotype given by 'observed_id' or just the first hap
    observed_hap = Haplotypes(observed_hap, log=log)
    observed_hap.read(haplotypes=(set((observed_id,)) if observed_id is not None else None))

    # load the genotypes
    variants = {v.id for v in causal_hap.variants}
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
    causal_gt = causal_hap.transform(gts).sum(axis=1)
    observed_gt = observed_hap.transform(gts)
    observed_ld = np.array([
        pearson_corr_ld(causal_gt, observed_gt.data[:, o].sum(axis=1))
        for o in range(observed_gt.data.shape[1])
    ])

    # return the strongest LD among all haps
    best_observed_hap_idx = np.abs(observed_ld).argmax()
    return observed_ld[best_observed_hap_idx]


def plot_params(params: npt.NDArray, vals: npt.NDArray, val_title: str, hap_counts: npt.NDArray = None):
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
    hap_counts: npt.NDArray, optional
        The number of haplotypes present in each hap file
    """
    also_plot_counts = hap_counts is not None
    fig, axs = plt.subplots(nrows=len(params.dtype)+1+also_plot_counts, ncols=1, sharex=True)
    # create a plot for the vals, first
    axs[0].plot(vals, "g-")
    axs[0].set_ylabel(val_title)
    axs[0].set_xticklabels([])
    if also_plot_counts:
        axs[1].plot(hap_counts, "m-")
        axs[1].set_ylabel("num_haps")
    # now, plot each of the parameter values on the other axes
    for idx, param in enumerate(params.dtype.names):
        axs[idx+1+also_plot_counts].plot(params[param], "-")
        axs[idx+1+also_plot_counts].set_ylabel(param)
    return fig


def plot_params_simple(params: npt.NDArray, vals: npt.NDArray, val_title: str, hap_counts: npt.NDArray = None):
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
    hap_counts: npt.NDArray, optional
        The number of haplotypes present in each hap file
    """
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True)
    # plot params against vals
    ax.plot(params, vals,  "g-")
    ax.set_xlabel(params.dtype.names[0])
    ax.set_ylabel(val_title, color="g")
    # also plot hap_counts as a second y-axis
    if hap_counts is not None:
        ax = ax.twinx()
        ax.plot(params, hap_counts, "m-")
        ax.set_ylabel("Number of observed haplotypes", color="m")
    return fig


@click.command()
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("observed_hap", type=click.Path(path_type=Path))
@click.argument("causal_hap", type=click.Path(path_type=Path))
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
    "-n",
    "--num-haps",
    is_flag=True,
    default=True,
    show_default=True,
    help="Whether to also depict the total number of haplotypes in the file",
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
    region: str = None,
    observed_id: str = None,
    causal_id: str = None,
    num_haps: bool = True,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "ERROR"
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

    # compute LD between the causal hap and the best observed hap across param vals
    ld_vals = np.abs(np.array([
        get_best_ld(
            gts,
            Path(str(observed_hap).format(**dict(zip(dtypes.keys(), params[idx])))),
            Path(str(causal_hap).format(**dict(zip(dtypes.keys(), params[idx])))),
            region=region,
            observed_id=observed_id,
            causal_id=causal_id,
            log=log
        )
        for idx in range(len(params))
    ]))

    hap_counts = None
    if num_haps:
        hap_counts = np.array([
            get_num_haps(Path(str(observed_hap).format(**dict(zip(dtypes.keys(), params[idx])))))
            for idx in range(len(params))
        ])

    if len(dtypes) > 1:
        fig = plot_params(params, ld_vals, "LD with causal hap", hap_counts)
    elif len(dtypes) == 1:
        fig = plot_params_simple(params, ld_vals, "LD with causal hap", hap_counts)
    else:
        raise ValueError("No parameter values found")

    fig.tight_layout()
    fig.savefig(output)

if __name__ == "__main__":
    main()
