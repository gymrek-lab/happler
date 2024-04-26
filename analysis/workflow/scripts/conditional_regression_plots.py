#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
import matplotlib
import numpy as np
matplotlib.use('Agg')
import numpy.typing as npt
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
    """
    pvals = -np.log10(pvals)
    if exclude_mask is not None:
        positions = positions[exclude_mask]
        pvals = pvals[exclude_mask]
        red_mask = red_mask[exclude_mask]
    if red_mask is not None:
        blue_mask = np.logical_not(red_mask)
        ax.scatter(positions[blue_mask], pvals[blue_mask])
        ax.scatter(positions[red_mask], pvals[red_mask], c="red")
    else:
        ax.scatter(positions, pvals)

def condition_on_variable(gts: Genotypes, pt: Phenotypes, covars: Genotypes = None):
    """
    Compute explained variance for each SNP or haplotype in a set

    The explained variance is beta^2 in a linear model y = bx + e when x and y have
    been standardized to mean 0 and stdev 1

    Parameters
    ----------
    pt: Phenotypes
        A Phenotypes object with only a single phenotype
    gts: Genotypes
        The genotypes of all of the SNPs
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
    return assoc_test.run(gts.data.sum(axis=2), pt.data[:, 0]).data["pval"]


@click.command()
@click.argument("genotypes", type=click.Path(path_type=Path))
@click.argument("phenotype", type=click.Path(path_type=Path))
@click.argument("haplotype", type=click.Path(path_type=Path))
@click.option(
    "-i",
    "--hap-id",
    type=str,
    show_default="the first haplotype",
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
    "--show-original",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to also depict the original Manhattan plot",
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
    show_original: bool = False,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Make Manhattan plots for a set of SNPs for each of three cases:
    1) the haplotype as a covariate
    2) the haplotypes' alleles as (multiple) covariates
    3) each of the haplotypes' alleles as a covariate, where we produce one plot for each allele
    """
    log = getLogger("plot_conditional_regressions", verbosity)

    # load a haplotype and its alleles
    hps = Haplotypes(haplotype, log=log)
    hps.read(haplotypes=(set((hap_id,)) if hap_id is not None else None))
    # if the user didn't specify a hap ID, just use the first one
    if hap_id is None:
        hap_id = next(iter(hps.data.keys()))
        hps.subset(haplotypes=(hap_id,), inplace=True)
    hap = hps.data[hap_id]
    # get the variants from the haplotype
    variants = tuple(var.id for var in hap.variants)

    # load just the first phenotype
    pts = Phenotypes(phenotype, log=log)
    pts.read()

    # load the SNPs
    gts = GenotypesVCF
    if genotypes.suffix == ".pgen":
        gts = GenotypesPLINK
    gts = gts(genotypes, log=log)
    gts.read(region=region, samples=set(pts.samples))
    # we need phasing for transform (below)
    gts.check_phase()
    gts.check_missing()
    gts.check_biallelic()
    positions = gts.variants["pos"]
    # reorder the phenotypes and ensure there is only one phenotype
    pts.subset(names=(pts.names[0],), samples=gts.samples, inplace=True)
    pts.standardize()

    # get the haplotype genotypes
    hap_gt = hps.transform(gts)

    # make the figure
    # set all panels in the same row
    figsize = (FIGSIZE*(len(variants)+2+show_original)/2.5, FIGSIZE)
    fig, axs = plt.subplots(1, 2+len(variants)+show_original, sharey=True, figsize=figsize)

    # highlight alleles in red
    red_mask = np.zeros(len(gts.variants), dtype=np.bool_)
    for snp in variants:
        red_mask[gts._var_idx[snp]] = True

    log.info("Creating haplotype plot")
    # first, encode the haplotype as covariate
    make_manhattan(axs[0], positions, condition_on_variable(gts, pts, hap_gt), red_mask)
    axs[0].set_title("Haplotype")
    log.info("Creating haplotype alleles plot")
    # now, encode the haplotypes' alleles as separate covariates
    covars = gts.subset(variants=variants)
    # exclude the covariates from plotting
    exclude = np.ones(len(gts.variants), dtype=np.bool_)
    for snp in variants:
        exclude[gts._var_idx[snp]] = False
    make_manhattan(axs[1], positions, condition_on_variable(gts, pts, covars), red_mask, exclude)
    axs[1].set_title("Haplotype's Alleles")
    # finally, encode each of the alleles as a covariate in a separate plot
    for idx in range(len(variants)):
        log.info(f"Creating plot #{idx+3}")
        covars = gts.subset(variants=(variants[idx],))
        exclude = np.ones(len(gts.variants), dtype=np.bool_)
        exclude[gts._var_idx[variants[idx]]] = False
        make_manhattan(axs[idx+2], positions, condition_on_variable(gts, pts, covars), red_mask, exclude)
        axs[idx+2].set_title(variants[idx])
    if show_original:
        log.info("Creating original manhattan plot")
        make_manhattan(axs[-1], positions, condition_on_variable(gts, pts), red_mask)
        axs[-1].set_title("")

    # now, tidy up and save the plot
    log.info("Writing out plot")
    fig.supylabel("-log10(pval)")
    fig.supxlabel("Chromosomal Position")
    fig.tight_layout()
    fig.savefig(output)

if __name__ == "__main__":
    main()
