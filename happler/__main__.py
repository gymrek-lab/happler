#!/usr/bin/env python

import click
from pathlib import Path

from . import Data, Tree


@click.group()
@click.version_option()
def main():
    """
    happler: A haplotype-based fine-mapping method

    Test for associations between a trait and haplotypes (ie sets of correlated SNPs rather than individual SNPs)
    """
    pass


@main.command()
@click.argument("genotypes", type=Path)
@click.argument("phenotypes", type=Path)
def run(genotypes: Path, phenotypes: Path):
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
    """
    gt = Data.Genotypes.load(genotypes)
    ph = Data.Phenotypes.load(phenotypes)
    tree = Tree.TreeBuilder(gt, ph)
