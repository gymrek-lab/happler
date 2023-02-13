#!/usr/bin/env python

import sys
from pathlib import Path

import click
import numpy as np

import haptools
from haptools.data import GenotypesPLINK


@click.command()
@click.argument("hps", type=click.Path(exists=True, path_type=Path))
@click.argument("samples_file", type=click.File("r"))
def main(hps: Path, samples_file: Path):
    """
    Create a column in the GT matrix from the provided haplotype pseudogenotypes

    Output the column to stdout

    \f
    Parameters
    ----------
    hps: Path
        The path to the pgen file containing the haplotype pseudogenotypes
    samples_file: Path
        A file containing just each sample ID on a different line of the file
    """
    gts = GenotypesPLINK.load(hps)

    with samples_file as samps_file:
        samples = samps_file.read().splitlines()

    gts.subset(samples=samples, inplace=True)

    assert len(gts.variants) == 1
    pos = gts.variants["pos"][0]

    hap_gts = gts.data[:, 0, :].sum(axis=1).astype(np.uint8)
    np.savetxt(sys.stdout.buffer, hap_gts, fmt='%i', delimiter='\n', header=f"{pos}:2", comments="")


if __name__ == "__main__":
    main()
