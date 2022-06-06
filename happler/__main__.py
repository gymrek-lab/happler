#!/usr/bin/env python

import sys
import click
import logging
from pathlib import Path
from typing import TextIO
from typing import Union, Tuple

from . import tree
from haptools import data


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
def main():
    """
    happler: A haplotype-based fine-mapping method

    Test for associations between a trait and haplotypes (ie sets of correlated SNPs) rather than individual SNPs
    """
    pass


@main.command(context_settings=CONTEXT_SETTINGS)
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
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
    "-s",
    "--sample",
    "samples",
    type=str,
    multiple=True,
    show_default="all samples",
    help=(
        "A list of the samples to subset from the genotypes file (ex: '-s sample1 -s"
        " sample2')"
    ),
)
@click.option(
    "-S",
    "--samples-file",
    type=click.File("r"),
    show_default="all samples",
    help=(
        "A single column txt file containing a list of the samples (one per line) to"
        " subset from the genotypes file"
    ),
)
@click.option(
    "--discard-multiallelic",
    is_flag=True,
    show_default="do not discard multi-allelic variants",
    help="Whether to discard multi-allelic variants or just complain about them."
)
@click.option(
    "-o",
    "--output",
    type=click.File("w"),
    default="-",
    show_default="stdout",
    help="A .hap file describing the extracted haplotypes.",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="ERROR",
    show_default="only errors",
    help="The level of verbosity desired",
)
def run(
    genotypes: Path,
    phenotypes: Path,
    region: str = None,
    samples: Tuple[str] = tuple(),
    samples_file: Path = None,
    discard_multiallelic: bool = False,
    output: TextIO = sys.stdout,
    verbosity: str = 'CRITICAL',
):
    """
    Use the tool to find trait-associated haplotypes

    GENOTYPES must be formatted as VCFs and

    PHENOTYPES must be a tab-separated file containing two columns: sample ID and
    phenotype value

    Ex: happler run tests/data/simple.vcf tests/data/simple.tsv > simple.hap

    \f
    Parameters
    ----------
    genotypes : Path
        The path to the genotypes in VCF format
    phenotypes : Path
        The path to the phenotypes in TSV format. There should be no header lines.
    region : str, optional
        See documentation for :py:meth:`~.data.Genotypes.read`
    sample : Tuple[str], optional
        See documentation for :py:meth:`~.data.Genotypes.read`
    samples_file : Path, optional
        A single column txt file containing a list of the samples (one per line) to
        subset from the genotypes file
    """
    log = logging.getLogger("run")
    logging.basicConfig(
        format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)",
        level=verbosity,
    )
    # handle samples
    if samples and samples_file:
        raise click.UsageError(
            "You may only use one of --sample or --samples-file but not both."
        )
    if samples_file:
        with samples_file as samps_file:
            samples = samps_file.read().splitlines()
    elif samples:
        # needs to be converted from tuple to list
        samples = list(samples)
    else:
        samples = None
    # load data
    log.info("Loading genotypes")
    if discard_multiallelic:
        gt = data.Genotypes(genotypes)
        gt.read(region=region, samples=samples)
        log.info("Discarding multiallelic variants")
        gt.check_biallelic(discard_also=True)
        gt.check_phase()
    else:
        gt = data.Genotypes.load(genotypes, region=region, samples=samples)
    log.info("There are {} samples and {} variants".format(*gt.data.shape))
    log.info("Loading phenotypes")
    ph = data.Phenotypes.load(phenotypes, samples=samples)
    log.info("Running tree builder")
    hap_tree = tree.TreeBuilder(gt, ph).run()
    log.info("Outputting haplotypes")
    tree.Haplotypes.from_tree(hap_tree).write(output)


if __name__ == "__main__":
    # run the CLI if someone tries 'python -m happler' on the command line
    main(prog_name="happler")
