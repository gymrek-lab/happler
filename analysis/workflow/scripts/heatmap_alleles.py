#!/usr/bin/env python
import click
import dataclasses
from pathlib import Path

# Allow us to edit fonts in Illustrator
import matplotlib
import numpy as np
matplotlib.use('Agg')
from haptools import logging
import matplotlib.pyplot as plt
from haptools.data import GenotypesVCF, GenotypesPLINK, Haplotypes, Phenotypes


def plot_hapmatrix(hpmt, hap_id, snps):
    """
    Adapted from this script by Shubham Saini
    https://github.com/shubhamsaini/pgc_analysis/blob/master/viz-snp-hap.py
    """
    box_w =  1.0/hpmt.shape[1]
    box_h = box_w
    hap_height = hpmt.shape[0]*0.0025*4
    legend_height = 0
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Plot SNPs
    ax.imshow(hpmt, vmin=0, vmax=1, cmap=plt.cm.Greys, aspect="auto", interpolation="none")
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xticks(np.arange(0, len(snps), 1))
    ax.set_xticklabels(snps)
    ax.set_title("All haplotypes" if hap_id == "ALL" else "Haplotype %s"%hap_id)
    return fig


@click.command()
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("phenotypes", type=click.Path(exists=True, path_type=Path))
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
        "A haplotype ID from the .hap file to plot"
        "(ex: '-i H1')."
    ),
)
@click.option(
    "--use-hap-alleles",
    is_flag=True,
    default=False,
    show_default=True,
    help="Color alleles from the haplotype as black instead of just the ALT alleles"
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
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def main(genotypes: Path, haplotypes: Path, phenotypes: Path, region: str = None, hap_id: str = None, use_hap_alleles: bool = False, output: Path = Path("/dev/stdout"), verbosity: str = "INFO"):
    """
    Create a heatmap plot that visualizes a haplotype from a .hap file
    """
    log = logging.getLogger("heatmap_alleles", level=verbosity)

    gts = GenotypesVCF
    if genotypes.suffix == ".pgen":
        gts = GenotypesPLINK

    gts = gts(genotypes, log=log)

    hps = Haplotypes(haplotypes, log=log)
    hps.read(haplotypes=(set((hap_id,)) if hap_id is not None else None))

    # get the variants from all haplotypes
    hps_vars = {var.id: var.allele for hap in hps.data for var in hps.data[hap].variants}
    if hap_id is None:
        hap_id = "ALL"
        if use_hap_alleles:
            log.warning("A hap ID wasn't given. Setting --use-hap-alleles to False")
            use_hap_alleles = False

    pts = Phenotypes(phenotypes, log=log)
    pts.read()

    gts.read(variants=set(hps_vars.keys()), region=region, samples=pts.samples)
    gts.check_phase()
    gts.check_missing()
    gts.check_biallelic()
    gts.subset(variants=list(hps_vars.keys()), inplace=True)
    pts.subset(samples=gts.samples)
    num_samps, num_snps, _ = gts.data.shape

    if use_hap_alleles:
        for var_idx, var in enumerate(gts.variants):
            # is it the REF or ALT allele?
            allele = var["alleles"][1] == hps_vars[var["id"]]
            # modify the genotype matrix to contain 1s only when the allele is present
            gts.data[:, var_idx] = gts.data[:, var_idx] == allele
            gts.variants[var_idx]["id"] += f":{int(allele)}"

    # resulting shape is num_samps * 2 by num_SNPs
    hpmt = gts.data.transpose((0, 2, 1)).reshape((num_samps*2, num_snps)).astype(np.bool_)
    # sort by each SNP from left to right
    samp_indices = np.lexsort(tuple(hpmt[:, i] for i in range(hpmt.shape[1]-1, -1, -1)))
    hpmt = hpmt[samp_indices]

    # also append the phenotypes
    pts = np.repeat(pts.data[:, 0], 2)[samp_indices]
    # standardize so that the max value is 1
    pts = (pts-pts.min())/(pts.max()-pts.min())
    hpmt = np.append(hpmt, pts[:, np.newaxis], axis=1)

    fig = plot_hapmatrix(hpmt, hap_id, list(gts.variants["id"])+["pheno"])
    fig.tight_layout()
    fig.savefig(output)


if __name__ == "__main__":
    main()
