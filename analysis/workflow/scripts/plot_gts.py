#!/usr/bin/env python
from pathlib import Path

import click
import numpy as np
import matplotlib.pyplot as plt

from haptools.logging import getLogger
from haptools.data import Phenotypes, GenotypesPLINK


@click.command()
@click.argument("gts", type=click.Path(exists=True, path_type=Path))
@click.argument("pts", type=click.Path(exists=True, path_type=Path))
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
    "--hap-id",
    type=str,
    show_default="all of the haplotypes",
    help=(
        "A haplotype ID from the gts file to plot"
        "(ex: '-i H1')."
    ),
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
    gts: Path,
    pts: Path,
    region: str = None,
    hap_id: str = None,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Plot phenotypes vs genotypes for a haplotype or multiple haplotypes
    """
    log = getLogger("plot_gts", verbosity)

    log.info("Loading haplotype genotypes")
    gts = GenotypesPLINK(fname=gts, log=log)
    gts.read(region=region, variants=set((hap_id,)) if hap_id is not None else None)
    gts.check_missing()
    gts.check_biallelic()
    hap_gts = gts.data[:, :, :2].sum(axis=2)

    pts = Phenotypes(pts, log=log)
    pts.read()
    pts.subset(samples=gts.samples)
    pt = pts.data[:, 0]

    # set up plot: each haplotype is a column in the subplots
    fig, axs = plt.subplots(1, len(gts.variants), sharey=True)

    for hp_idx in range(len(gts.variants)):
        # gt vals (should just be 0, 1, 2)
        gt_vals = np.unique(hap_gts[:, hp_idx])
        # group phenotype data by genotype
        grouped_phenos = []
        for gt in gt_vals:
            grouped_phenos.append(pt[hap_gts[:, hp_idx] == gt])

        # Plot box and whisker plot for each genotype value
        axs[hp_idx].boxplot(grouped_phenos, labels=gt_vals)
        axs[hp_idx].set_xlabel(gts.variants["id"][hp_idx])

    fig.supylabel("Phenotype")
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig(output)


if __name__ == "__main__":
    main()
