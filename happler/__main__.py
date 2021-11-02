#!/usr/bin/env python

import click


@click.group()
@click.version_option()
def main():
    """
    happler: A haplotype-based fine-mapping method\n
    Test for associations between a trait and haplotypes (sets of correlated variants) rather than individual SNPs
    """
    pass


@main.command()
def run():
    """
    Use the tool to find trait-associated haplotypes
    """
    pass
