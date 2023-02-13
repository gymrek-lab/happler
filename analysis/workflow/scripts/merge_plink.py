#!/usr/bin/env python

import click
import numpy as np
from pathlib import Path
from haptools.data import GenotypesPLINK


@click.command()
@click.argument("file1", type=click.Path(exists=True, path_type=Path))
@click.argument("file2", type=click.Path(exists=True, path_type=Path))
@click.argument("output", type=click.Path(path_type=Path))
@click.option(
    "--replace/--no-replace",
    is_flag=True,
    default=True,
    show_default=True,
    help="Whether to use the variants in file2 to replace those in file1",
)
def main(file1: Path, file2: Path, output: Path, replace: bool = True):
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
    replace: bool, optional
        If True, take the left-inner join of the set of variants, so that
        conflicting IDs from file2 are used to replace the GTs in file1

        You should only set this flag to False if you are confident that file1 and
        file2 have no conflicting variant IDs. Otherwise, you may end up with an
        output fileset that has duplicate IDs!
    """
    gts1 = GenotypesPLINK(file1)
    gts2 = GenotypesPLINK(file2)

    gts1.read()
    gts2.read()

    gts1.check_phase()
    gts2.check_phase()

    assert gts1.samples == gts2.samples

    if replace:
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

    # append any variants from file2 to the end of file1
    gts1.variants = np.concatenate((gts1.variants, gts2.variants))
    gts1.data = np.concatenate((gts1.data, gts2.data), axis=1)

    gts1.fname = output
    gts1.write()
    

if __name__ == "__main__":
    main()
