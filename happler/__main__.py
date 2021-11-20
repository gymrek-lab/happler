#!/usr/bin/env python

import click
from pathlib import Path
from typing import Union, Tuple

from . import data, tree


@click.group()
@click.version_option()
def main():
    """
    happler: A haplotype-based fine-mapping method

    Test for associations between a trait and haplotypes (ie sets of correlated SNPs rather than individual SNPs)
    """
    pass


@main.command()
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
def run(
    genotypes: Path,
    phenotypes: Path,
    region: str = None,
    samples: Tuple[str] = tuple(),
    samples_file: Path = None,
):
    """
    Use the tool to find trait-associated haplotypes

    GENOTYPES must be formatted as VCFs and
    PHENOTYPES must be a tab-separated file containing two columns: sample ID and
    phenotype value

    \f
    Parameters
    ----------
    genotypes : Path
        The path to the genotypes in VCF format
    phenotypes : Path
        The path to the phenotypes in TSV format
    region : str, optional
        See documentation for :py:meth:`~.data.Genotypes.read`
    sample : Tuple[str], optional
        See documentation for :py:meth:`~.data.Genotypes.read`
    samples_file : Path, optional
        A single column txt file containing a list of the samples (one per line) to
        subset from the genotypes file
    """
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
    gt = data.Genotypes.load(genotypes, region=region, samples=samples)
    ph = data.Phenotypes.load(phenotypes, samples=samples)
    hap_tree = tree.TreeBuilder(gt, ph)

if __name__ == '__main__':
    # run the CLI if someone tries 'python -m happler' on the command line
    main(prog_name="happler")
