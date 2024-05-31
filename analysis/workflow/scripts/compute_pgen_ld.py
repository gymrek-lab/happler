#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
import numpy as np

from haptools.logging import getLogger
from haptools.ld import pearson_corr_ld
from haptools.data import Data, Genotypes, GenotypesPLINK


def corr(a, b):
    """
    Speedily compute the pearson correlation coefficient of a set of variables
    against a single variable

    Adapted from https://stackoverflow.com/a/61587017/16815703

    Please note that I haven't checked whether this offers the same precision as
    np.corrcoef, but it is certainly faster.

    Parameters
    ----------
        a: 2d numpy array
        b: 1d numpy array

    Returns
    -------
        c: numpy array
            correlation coefficients of all columns of a against b
    """
    a = np.append(a, b[:, np.newaxis], axis=1)
    i = -1

    mean_t = np.mean(a, axis=0)
    std_t = np.std(a, axis=0)
    mean_i = mean_t[i]
    std_i = std_t[i]
    mean_xy = np.mean(a*a[:,i][:,None], axis=0)
    c = (mean_xy - mean_i * mean_t)/(std_i * std_t)
    return c[:-1]


@click.command()
@click.argument("gts", type=click.Path(exists=True, path_type=Path))
@click.argument("target", type=click.Path(exists=True, path_type=Path))
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
    show_default="all SVs",
    help="Only select SVs with a MAF above this threshold",
)
@click.option(
    "-i",
    "--hap-id",
    type=str,
    show_default="the first haplotype ID",
    help=(
        "A haplotype ID from the target file to use "
        "(ex: '-i H1')."
    ),
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A .ld file containing the LD results for each SNP in gts",
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
    target: Path,
    region: str = None,
    maf: float = None,
    hap_id: str = None,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Compute LD between a SNP in a PGEN file and all other SNPs in a different
    PGEN file
    """
    log = getLogger("compute_pgen_ld", verbosity)

    log.info("Loading target genotypes")
    target = GenotypesPLINK(fname=target, log=log)
    target.read(variants=set((hap_id,)) if hap_id is not None else None)
    target.check_missing()
    target.check_biallelic()

    if hap_id is None:
        target.subset(variants=(target.variants["id"][0],), inplace=True)

    log.info("Loading reference genotypes")
    gts = GenotypesPLINK(fname=gts, log=log)
    gts.read(samples=set(target.samples), region=region)
    gts.check_missing()
    gts.check_biallelic()
    gts.check_maf(threshold=maf, discard_also=True)
    # important: check that samples are ordered the same in each file!
    assert gts.samples == target.samples

    log.info("Summing target genotypes")
    target_gts = target.data[:, 0, :2].sum(axis=1)
    variant_gts = gts.data[:, :, :2].sum(axis=2)

    variant_lds = corr(variant_gts, target_gts)

    log.info("Computing LD between genotypes and the target")
    with Data.hook_compressed(output, mode="w") as ld_file:
        log.info("Outputting .ld file with LD values")
        ld_file.write("CHR\tBP\tSNP\tR\n")
        for idx, variant in enumerate(gts.variants):
            var_chr, var_bp, var_snp = variant[["chrom", "pos", "id"]]
            variant_ld = variant_lds[idx]
            ld_file.write(f"{var_chr}\t{var_bp}\t{var_snp}\t{variant_ld:.3f}\n")


if __name__ == "__main__":
    main()
