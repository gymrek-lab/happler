#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
from haptools.logging import getLogger
from haptools.sim_phenotype import Haplotype
from haptools.data import GenotypesPLINK, Phenotypes, Haplotypes

from happler.tree.assoc_test import AssocTestSimpleSMTScore

@click.command()
@click.argument("gt", type=click.Path(exists=True, path_type=Path))
@click.argument("pt", type=click.Path(exists=True, path_type=Path))
@click.argument("hp", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A transformed genotypes file",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def main(
    gt: Path,
    pt: Path,
    hp: Path,
    output: Path,
    verbosity: str = "DEBUG",
):
    """
    Rewrite a two-allele haplotype in a .hap file. Output the new .hap file to stdout.

    If the BIC of the second allele is less than the first, swap them so that the second allele becomes first.

    Otherwise, just don't print anything.
    """
    log = getLogger("flip_hap_alleles", verbosity)

    hap = Haplotypes(hp, haplotype=Haplotype)
    hap.read()

    hap_vars = tuple(v.id for v in list(hap.data.values())[0].variants)
    phen = Phenotypes(pt)
    phen.read()

    hap_gts = GenotypesPLINK.load(gt)
    hap_gts = hap_gts.subset(variants=hap_vars)

    data = AssocTestSimpleSMTScore(with_bic=True).run(
        hap_gts.data.sum(axis=2),
        phen.data[:, 0],
    ).data

    if data["bic"][1] < data["bic"][0]:
        hap.data["H0"].variants = hap.data["H0"].variants[::-1]
        hap.fname = output
        hap.write()


if __name__ == "__main__":
    main()
