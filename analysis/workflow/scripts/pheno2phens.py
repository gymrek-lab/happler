#!/usr/bin/env python

import click
import numpy as np
from pathlib import Path

import haptools
from haptools.data import Data, GenotypesPLINK, Phenotypes


@click.command()
@click.argument("hps", type=click.Path(exists=True, path_type=Path))
@click.argument("pheno", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A TSV file containing simulated phenotypes",
)
def main(hps: Path, pheno: Path, output: Path):
    """
    Create a phens file for a causal haplotype

    \f
    Parameters
    ----------
    hps: Path
        The path to the pgen file containing the haplotype pseudogenotypes
    pheno: Path
        A file containing the pheno file for the causal haplotype
    output: Path
        The path to the output phens file
    """
    gts = GenotypesPLINK.load(hps)
    pts = Phenotypes.load(pheno)

    assert len(gts.variants) == 1
    assert gts.samples == pts.samples
    assert len(pts.names) == 1
    pos = gts.variants["pos"][0]
    
    hap_gts = gts.data[:, 0, :].sum(axis=1)
    # normalize the genotypes
    hap_gts = (hap_gts - hap_gts.mean()) / hap_gts.std()

    with Data.hook_compressed(output, mode="w") as phens:
        phens.write(f"sample\t{pos}:2\tphen\n")
        for samp, gt, pt in zip(gts.samples, hap_gts, pts.data[:, 0]):
            phens.write(f"{samp}\t{gt:.17f}\t{pt:.17f}\n")


if __name__ == "__main__":
    main()
