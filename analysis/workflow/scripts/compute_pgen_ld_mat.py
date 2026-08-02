#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
import numpy as np

from haptools.logging import getLogger
from haptools.ld import pearson_corr_ld
from haptools.data import Data, GenotypesPLINK, GenotypesPLINKTR


@click.command()
@click.argument("snps", type=click.Path(exists=True, path_type=Path))
@click.argument("repeats", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all genotypes",
    help="""
    The region from which to extract genotypes; ex: 'chr1:1234-34566' or 'chr7'\n
    For this to work, the seqnames must match!""",
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help=(
        "Perform reading operations in chunks of X variants. "
        "This reduces memory but at the cost of time."
    ),
)
@click.option(
    "--vcftype",
    type=str,
    default=None,
    show_default="infer",
    help="The type of repeat file. See https://haptools.readthedocs.io/en/stable/formats/genotypes.html#tandem-repeats",
)
@click.option(
    "--discard-missing",
    is_flag=True,
    show_default=True,
    default=False,
    help="Discard any variants with missing values instead of raising an error",
)
@click.option(
    "--maf",
    type=float,
    default=None,
    show_default="all SVs",
    help="Only select SNPs with a MAF above this threshold",
)
@click.option(
    "--r2",
    is_flag=True,
    show_default=True,
    default=False,
    help="Compute r^2 (unsigned) instead of r (signed)",
)
@click.option(
    "--precision",
    type=int,
    default=6,
    show_default=True,
    help="How many points after the decimal?",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A .vcor2 file for the LD matrix containing pairwise LD between all SNPs and repeats",
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
    snps: Path,
    repeats: Path,
    region: str = None,
    chunk_size: int = None,
    vcftype: str = None,
    discard_missing: bool = False,
    maf: float = None,
    r2: bool = False,
    precision: int = 6,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Compute pairwise LD between all SNPs and repeats in two different PGEN files

    Note that this command computes r (signed) not r^2 (unsigned). This corresponds w/
    plink2 --r-unphased 'square0' 'inter-chr' 'ref-based' 'yes-really' --nonfounders --ld-window-r2 0
    You can use --r2 to change this to match --r2-unphased
    """
    log = getLogger("compute_pgen_ld_mat", verbosity)

    log.info("Loading repeat genotypes")
    repeats = GenotypesPLINKTR(fname=repeats, chunk_size=chunk_size, log=log, vcftype=(vcftype if vcftype is not None else 'auto'))
    repeats.read(region=region)
    repeats.check_missing(discard_also=discard_missing)

    log.info("Loading SNP genotypes")
    snps = GenotypesPLINK(fname=snps, chunk_size=chunk_size, log=log)
    snps.read(samples=set(repeats.samples), region=region)
    snps.check_missing(discard_also=discard_missing)
    snps.check_biallelic()
    snps.check_maf(threshold=maf, discard_also=True)

    # important: check that samples are ordered the same in each file!
    repeats.subset(samples=snps.samples)
    assert snps.samples == repeats.samples

    all_var_ids = np.hstack((snps.variants["id"], repeats.variants["id"]))
    if output != Path("/dev/stdout"):
        log.info("Outputting .vars file listing variant IDs")
        var_ids_out = output + ".vars"
        if var_ids_out.suffix == ".vcor2.vars":
            var_ids_out = var_ids_out.with_suffix("").with_suffix(".vars")
        with Data.hook_compressed(var_ids_out, mode="w") as vars_file:
            vars_file.write("\n".join(all_var_ids))

    log.info("Summing and concatenating genotypes")
    repeats_gts = repeats.data[:, :, :2].sum(axis=2)
    snps_gts = snps.data[:, :, :2].sum(axis=2)
    all_gts = np.concatenate((snps_gts, repeats_gts), axis=1)

    log.info("Computing LD matrix")
    all_ld = np.corrcoef(all_gts, rowvar=False)

    log.info("Zeroing upper-right triangle in-place")
    upper_indices = np.triu_indices_from(all_ld, k=1)
    all_ld[upper_indices] = 0

    if r2:
        log.info("Squaring all values in-place")
        np.square(all_ld, out=all_ld)

    log.info("Outputting .vcor2 file with LD matrix")
    with Data.hook_compressed(output, mode="w") as ld_file:
        header = ''
        if output == Path("/dev/stdout"):
            header = "\t".join(all_var_ids)
        np.savetxt(ld_file, all_ld, delimiter='\t', fmt=f"%.{precision}f", header=header, comments='')


if __name__ == "__main__":
    main()
