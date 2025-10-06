#!/usr/bin/env python
import re
from pathlib import Path
from logging import Logger

import click
import matplotlib
import numpy as np
matplotlib.use('Agg')
import numpy.typing as npt
import statsmodels.api as sm
import matplotlib.pyplot as plt

from haptools.logging import getLogger
from happler.tree.assoc_test import AssocTestSimpleSM, AssocTestSimpleCovariates
from haptools.data import (
    Genotypes,
    Phenotypes,
    Haplotypes,
    GenotypesVCF,
    GenotypesPLINK,
)

FIGSIZE = 6


def make_manhattan(
    ax: plt.Axes,
    positions: npt.NDArray,
    pvals: npt.NDArray,
    red_mask: npt.NDArray = None,
    exclude_mask: npt.NDArray = None,
    orange_mask: npt.NDArray = None,
):
    """
    Add a Manhattan plot to the desired Axes object

    Parameters
    ----------
    ax: plt.Axes
        The matplotlib axes to plot onto
    positions: npt.NDArray
        A 1D array of base pair positions for each SNP
    pvals: npt.NDArray
        A 1D array of p-values of the same length as positions
    red_mask: npt.NDArray, optional
        If provided, a bool array of pvals to highlight in red
    exclude_mask: npt.NDArray, optional
        If provided, a bool array of pvals to exclude from plotting
    orange_mask: npt.NDArray, optional
    """
    try:
        pvals = -np.log10(pvals)
    except TypeError:
        # handle Decimals mixed in with np.float64s
        for idx in range(len(pvals)):
            try:
                pvals[idx] = -np.log10(pvals[idx])
            except TypeError:
                pvals[idx] = -np.float64(pvals[idx].log10())
    pvals = pvals.astype(np.float64)
    if exclude_mask is not None:
        positions = positions[exclude_mask]
        pvals = pvals[exclude_mask]
        red_mask = red_mask[exclude_mask]
    if red_mask is not None:
        blue_mask = np.logical_not(red_mask)
        ax.scatter(positions[blue_mask], pvals[blue_mask])
        ax.scatter(positions[red_mask], pvals[red_mask], c="red")
        if orange_mask is not None:
            ax.scatter(positions[orange_mask], pvals[orange_mask], c="orange")
    else:
        ax.scatter(positions, pvals)

def regress_effect(pt: Phenotypes, covars: Genotypes = None):
    """
    Regress out the effect(s) of the covariates on the phenotype

    Parameters
    ----------
    pt: Phenotypes
        A Phenotypes object with only a single phenotype
    covars: Genotypes, optional
        The variables on which to condition the SNPs by encoding them as covariates

    Returns
    -------
    Phenotypes
        The new phenotype values
    """
    resids = Phenotypes(fname=None, log=pt.log)
    resids.samples = pt.samples
    resids.names = pt.names
    if covars is None or not covars.data.shape[1]:
        resids.data = pt.data.copy()
    else:
        resids.data = (
            sm.OLS(pt.data, sm.add_constant(covars.data.sum(axis=2)))
            .fit()
            .resid[:, np.newaxis]
        )
    return resids

def condition_on_variable(gts: npt.NDArray, pt: Phenotypes, covars: Genotypes = None):
    """
    Compute pvals for each SNP or haplotype in a set after conditioning on others

    Parameters
    ----------
    gts: npt.NDArray
        The genotypes of all of the SNPs
    pt: Phenotypes
        A Phenotypes object with only a single phenotype
    covars: Genotypes, optional
        The variables on which to condition the SNPs by encoding them as covariates

    Returns
    -------
    npt.NDArray[float]
        The p-values of all of the input SNPs when conditioned on the covariates
    """
    if covars is None:
        assoc_test = AssocTestSimpleSM()
    else:
        assoc_test = AssocTestSimpleCovariates(covars=covars.data.sum(axis=2))
    return assoc_test.run(gts.sum(axis=2), pt.data[:, 0]).data["pval"]

def condition_on_variable_chunked(
    gts: Genotypes,
    pt: Phenotypes,
    covars: Genotypes = None,
    chunk_size: int = None,
):
    pvals = np.empty(gts.data.shape[1], dtype=object)

    chunks = chunk_size
    if chunks is None or chunks > len(pvals):
        chunks = len(pvals)

    for start in range(0, len(pvals), chunks):
        end = start + chunks
        if end > len(pvals):
            end = len(pvals)
        size = end - start

        pvals[start:end] = condition_on_variable(gts.data[:, start:end], pt, covars)

    return pvals

def extract_snp_matrix(logfile: Path) -> npt.NDArray:

    text = logfile.read_text().splitlines()
    iterations = {}
    current_iter = None
    current_tree = None

    for lineno, line in enumerate(text, start=1):
        ls = line.strip()
        # robustly find "Iteration X: tree Y" even if line is indented
        m = re.search(r"Iteration\s+(\d+)\s*:\s*tree\s+(\d+)", ls, re.I)
        if m:
            current_iter, current_tree = map(int, m.groups())
            continue

        # find label lines (skip lines mentioning 'root')
        if '[label="' in line and 'root' not in line:
            # try to capture up to the literal backslash-n ("\\n") first
            m2 = re.search(r'\[label="([^\\]+?)\\n', line)
            if not m2:
                # fallback: capture up to the next closing quote
                m2 = re.search(r'\[label="([^"]+)"', line)
            if not m2:
                raise ValueError(f"Couldn't parse [label=...] on line {lineno}: {line!r}")
            snp = m2.group(1)

            if current_iter is None or current_tree is None:
                raise ValueError(
                    f"Found SNP label on line {lineno} but no preceding 'Iteration ...: tree ...' line."
                )

            # if there are multiple SNPs in a tree, just choose the first one
            if current_iter not in iterations or current_tree not in iterations[current_iter]:
                iterations.setdefault(current_iter, {})[current_tree] = snp

    if not iterations:
        return np.empty((0, 0), dtype=object)

    # enforce same tree keys for every iteration
    first_it = sorted(iterations.keys())[0]
    expected_keys = sorted(iterations[first_it].keys())
    for it in sorted(iterations.keys()):
        keys = sorted(iterations[it].keys())
        if keys != expected_keys:
            raise ValueError(f"Inconsistent tree keys for iteration {it}: {keys} != {expected_keys}")

    rows = []
    for it in sorted(iterations.keys()):
        row = [iterations[it][k] for k in expected_keys]
        rows.append(row)

    return np.array(rows, dtype=object)


@click.command()
@click.argument("genotypes", type=click.Path(path_type=Path))
@click.argument("phenotype", type=click.Path(path_type=Path))
@click.argument("haplotype", type=click.Path(path_type=Path))
@click.option(
    "-i",
    "--hap-id",
    type=str,
    show_default="all haplotypes",
    help=(
        "A haplotype ID from the .hap file to plot"
        "(ex: '-i H1')."
    ),
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
    "--maf",
    type=float,
    default=None,
    show_default="all SNPs",
    help="Only select SNPs with a MAF above this threshold",
)
@click.option(
    "--show-original",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to also depict the original Manhattan plot",
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help=(
        "Perform reading/writing operations in chunks of X variants. "
        "This reduces memory but at the cost of time."
    ),
)
@click.option(
    "--log-file",
    type=Path,
    default=None,
    show_default="no extra plots",
    help=(
        "If provided, we will read trees from this happler log file and attempt to "
        "plot them as additional rows in the conditional regression plot"
    )
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A .png file containing the output plot",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="DEBUG",
    show_default=True,
    help="The level of verbosity desired",
)
def main(
    genotypes: Path,
    phenotype: Path,
    haplotype: Path,
    hap_id: str = None,
    region: str = None,
    maf: float = None,
    show_original: bool = False,
    chunk_size: int = None,
    log_file: Path = None,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Make Manhattan plots for a set of SNPs for each of three cases:
    1) the haplotype as a covariate
    2) the haplotypes' alleles as (multiple) covariates
    3) each of the haplotypes' alleles as a covariate, where we produce one plot for each allele
    4) the original manhattan plot without any conditioning
    """
    log = getLogger("plot_conditional_regressions", verbosity)

    # load a haplotype and its alleles
    hps = Haplotypes(haplotype, log=log)
    hps.read(haplotypes=(set((hap_id,)) if hap_id is not None else None))
    # get the variants from the haplotype
    # use dict to remove duplicate IDs
    variants = tuple(dict.fromkeys(var.id for hap in hps.data.values() for var in hap.variants))

    # load just the first phenotype
    pts = Phenotypes(phenotype, log=log)
    pts.read()

    # load the SNPs
    gts = GenotypesVCF
    if genotypes.suffix == ".pgen":
        gts = GenotypesPLINK
    gts = gts(genotypes, log=log, chunk_size=chunk_size)
    gts.read(region=region, samples=set(pts.samples))
    # we need phasing for transform (below)
    gts.check_phase()
    gts.check_missing()
    gts.check_biallelic()
    gts.check_maf(threshold=maf, discard_also=True)
    positions = gts.variants["pos"]
    # reorder the phenotypes and ensure there is only one phenotype
    pts.subset(names=(pts.names[0],), samples=gts.samples, inplace=True)
    pts.standardize()

    # get the haplotype genotypes
    hap_gt = hps.transform(gts)

    # read SNPs from the log file
    other_alleles = np.array([], dtype=object)
    if log_file is None:
        # make the figure
        # set all panels in the same row
        figsize = (FIGSIZE*(len(variants)+2+show_original)/2.5, FIGSIZE)
        fig, axs = plt.subplots(1, 2+len(variants)+show_original, sharey=True, figsize=figsize)
    else:
        other_alleles = extract_snp_matrix(log_file)

        # make the figure
        # set all panels in the same row
        figsize = (FIGSIZE*(len(variants)+2+show_original)/2.5, FIGSIZE*other_alleles.shape[0])
        fig, axs = plt.subplots(1+other_alleles.shape[0], 2+len(variants)+show_original, sharey=True, figsize=figsize)
        all_axs = axs
        axs = axs[0]

    # highlight alleles in red
    red_mask = np.zeros(len(gts.variants), dtype=np.bool_)
    for snp in variants:
        red_mask[gts._var_idx[snp]] = True

    log.info("Creating haplotype plot")
    # first, encode the haplotype as covariate
    make_manhattan(
        axs[0],
        positions,
        condition_on_variable_chunked(gts, pts, hap_gt, chunk_size=chunk_size),
        red_mask
    )
    axs[0].set_title("Haplotype"+("s" if len(hap_gt.variants)-1 else ""))
    log.info("Creating haplotype alleles plot")
    # now, encode the haplotypes' alleles as separate covariates
    covars = gts.subset(variants=variants)
    # exclude the covariates from plotting
    exclude = np.ones(len(gts.variants), dtype=np.bool_)
    for snp in variants:
        exclude[gts._var_idx[snp]] = False
    make_manhattan(
        axs[1],
        positions,
        condition_on_variable_chunked(gts, pts, covars, chunk_size=chunk_size),
        red_mask,
        exclude,
    )
    # f-string expressions cannot include a backslash
    f_str_quotation_plural_agh = 's\'' if len(hap_gt.variants)-1 else '\'s'
    axs[1].set_title(f"Haplotype{f_str_quotation_plural_agh} Alleles")
    # finally, encode each of the alleles as a covariate in a separate plot
    for idx in range(len(variants)):
        log.info(f"Creating plot #{idx+3}")
        covars = gts.subset(variants=(variants[idx],))
        exclude = np.ones(len(gts.variants), dtype=np.bool_)
        exclude[gts._var_idx[variants[idx]]] = False
        make_manhattan(
            axs[idx+2],
            positions,
            condition_on_variable_chunked(gts, pts, covars, chunk_size=chunk_size),
            red_mask,
            exclude,
        )
        axs[idx+2].set_title(variants[idx])
    if show_original:
        log.info("Creating original manhattan plot")
        # we're going to have to append the haplotypes to the original SNP genotypes,
        # so we'll have to extend the red mask to include the haplotypes
        # and create an orange mask for the haplotypes themselves
        red_mask = np.append(red_mask, np.zeros(len(hap_gt.variants), dtype=np.bool_))
        orange_mask = np.zeros(red_mask.shape[0], dtype=np.bool_)
        orange_mask[-len(hap_gt.variants):] = 1
        new_gt = Genotypes.merge_variants((gts, hap_gt), fname=None, log=log)
        # now, let's finally plot everything and append the haps to the gts
        make_manhattan(
            axs[-1],
            new_gt.variants["pos"],
            condition_on_variable_chunked(new_gt, pts, chunk_size=chunk_size),
            red_mask,
            exclude_mask=None,
            orange_mask=orange_mask,
        )
        axs[-1].set_title("")
    
    if log_file is not None:
        for iteration_idx in range(other_alleles.shape[0]):
            for tree_idx in range(other_alleles.shape[1]):
                if iteration_idx:
                    # regress out all but the current tree from the last iteration
                    covar_variants = other_alleles[iteration_idx-1].tolist()
                    covar_variants.remove(other_alleles[iteration_idx-1, tree_idx])
                else:
                    # if this is the first iteration, we need to do special things
                    covar_variants = other_alleles[0, :tree_idx].tolist()
                # highlight the variant that we found in red
                red_mask = np.zeros(len(gts.variants), dtype=np.bool_)
                red_mask[gts._var_idx[other_alleles[iteration_idx, tree_idx]]] = True
                # set title to indicate the variants that we are conditioning on
                all_axs[iteration_idx+1, tree_idx].set_title("\n".join(covar_variants))
                covars = gts.subset(variants=set(covar_variants))
                resids = regress_effect(pts, covars)
                make_manhattan(
                    all_axs[iteration_idx+1, tree_idx],
                    positions,
                    condition_on_variable_chunked(gts, resids, chunk_size=chunk_size),
                    exclude_mask=None,
                    red_mask=red_mask,
                    orange_mask=None,
                )
            for axes_idx in range(other_alleles.shape[1], all_axs.shape[1]):
                all_axs[iteration_idx+1, axes_idx].remove()

    # now, tidy up and save the plot
    log.info("Writing out plot")
    fig.supylabel("-log10(pval)")
    fig.supxlabel("Chromosomal Position")
    fig.suptitle("Manhattan plots when conditioning on...")
    fig.tight_layout()
    fig.savefig(output)

if __name__ == "__main__":
    main()
