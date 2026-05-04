#!/usr/bin/env python
from pathlib import Path

import click
import numpy as np

from haptools.logging import getLogger
from haptools.data import GenotypesPLINK


@click.command()
@click.argument("file1", type=click.Path(exists=True, path_type=Path))
@click.argument("file2", type=click.Path(exists=True, path_type=Path))
@click.argument("output", type=click.Path(path_type=Path))
@click.option(
    "--region",
    default=None,
    show_default=True,
    help="The region to extract genotypes from",
)
@click.option(
    "--maf",
    type=float,
    default=None,
    show_default="all variants",
    help="Only use variants with MAFs above this threshold for both files",
)
@click.option(
    "--maf-file",
    type=click.Choice(["1", "2", "both"]),
    default="both",
    show_default=True,
    help="Which file should the MAF threshold be applied to?",
)
@click.option(
    "--replace/--no-replace",
    is_flag=True,
    default=True,
    show_default=True,
    help="Whether to use the variants in file2 to replace those in file1",
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
    "--extract",
    type=click.File("r"),
    default=None,
    show_default="all variants",
    help="Keep only certain variants from file1. Works similarly to plink2's --extract",
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
    file1: Path,
    file2: Path,
    output: Path,
    region: str = None,
    maf: float = None,
    maf_file: str = "both",
    replace: bool = True,
    chunk_size: int = None,
    extract: Path = None,
    verbosity: str = "DEBUG",
):
    """
    Merge variants from two PGEN files that have the same set of samples

    \f
    Parameters
    ----------
    file1: Path
        The path to the first pgen file
    file2: Path
        The path to the second pgen file
    output: Path
        The path to the output pgen file
    region: Path
        The region to extract genotypes from
    replace: bool, optional
        If True, take the left-inner join of the set of variants, so that
        conflicting IDs from file2 are used to replace the GTs in file1

        You should only set this flag to False if you are confident that file1 and
        file2 have no conflicting variant IDs. Otherwise, you may end up with an
        output fileset that has duplicate IDs!
    extract: Path, optional
        This will do the same thing as plink2's --extract argument but for file1 only
    verbosity: str, optional
        How verbose do we want the log to be?
    """
    log = getLogger("merge_plink", verbosity)

    log.info("Loading genotypes from both files")
    gts1 = GenotypesPLINK(fname=file1, chunk_size=chunk_size, log=log)
    gts2 = GenotypesPLINK(fname=file2, chunk_size=chunk_size, log=log)

    if extract is not None:
        with extract as extract_file:
            extract = set(extract_file.read().splitlines())

    gts1.read(region=region, variants=extract)
    gts2.read(region=region)

    if maf is not None:
        if maf_file == "both":
            maf_file = "both files"
        elif maf_file == "1":
            maf_file = "file1"
        else:
            maf_file = "file2"
        log.info(f"Subsetting {maf_file} by MAF")
        if maf_file in ("file1", "both files"):
            gts1.check_missing()
            gts1.check_biallelic()
            gts1.check_maf(threshold=maf, discard_also=True)
        if maf_file in ("file2", "both files"):
            gts2.check_missing()
            gts2.check_biallelic()
            gts2.check_maf(threshold=maf, discard_also=True)

    if gts1.samples != gts2.samples:
        log.info("Getting intersection of samples in order of file1")
        samples = frozenset(gts2.samples)
        samples = tuple(samp for samp in gts1.samples if samp in samples)
        for gts in (gts1, gts2):
            gts.subset(samples=samples, inplace=True)

    if replace:
        log.info("Replacing variants with shared IDs")
        # which variants are shared? what are their indices within each file?
        common_ids, idxs_in_1, idxs_in_2 = np.intersect1d(
            gts1.variants["id"],
            gts2.variants["id"],
            return_indices=True
        )
        # update variants in file1 with their counterparts in file2
        gts1.variants[idxs_in_1] = gts2.variants[idxs_in_2]
        gts1.data[:, idxs_in_1, :] = gts2.data[:, idxs_in_2, :]
        # remove the common variants from file2
        gts2.variants = np.delete(gts2.variants, idxs_in_2)
        gts2.data = np.delete(gts2.data, idxs_in_2, axis=1)

    log.info("Appending any variants from file2 to the end of file1")
    gts1.variants = np.concatenate((gts1.variants, gts2.variants))
    gts1.data = np.concatenate((gts1.data, gts2.data), axis=1)

    log.info("Writing output")
    gts1.fname = output
    gts1.write()
    

if __name__ == "__main__":
    main()
