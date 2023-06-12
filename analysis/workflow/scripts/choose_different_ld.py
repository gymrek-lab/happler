#!/usr/bin/env python

import sys
from pathlib import Path
from logging import Logger
from itertools import combinations, product

import click
import numpy as np

from haptools.logging import getLogger
from haptools.ld import pearson_corr_ld
from haptools.simphenotype import Haplotype
from haptools.data import Genotypes, GenotypesPLINK, Haplotypes, Variant


def find_haps(gts: Genotypes, log: Logger, min_ld: float, max_ld: float, reps: float = 1, step: float = 0.1):
    log.info(f"Using min_ld {min_ld}, max_ld {max_ld}, reps {reps}, and step {step}")
    sum_gts = gts.data.sum(axis=2)

    chrom = np.unique(gts.variants["chrom"])
    assert len(chrom) == 1
    chrom = chrom[0]

    intervals = np.arange(max(0, min_ld), min(1, max_ld), step) + step
    ld_bins = {idx: 0 for idx in range(len(intervals))}
    count = len(intervals) * reps

    combos = np.array(list(range(len(gts.variants))))
    np.random.shuffle(combos)
    combo_ids = combinations(combos, 2)
    for combo_id in combo_ids:
        snp_vars = gts.variants[list(combo_id)]
        hp_id = "-".join(snp_vars['id'])
        log.debug(f"Trying {hp_id}")
        pos1, pos2 = snp_vars["pos"]
        for allele_combo in product(*snp_vars["alleles"]):
            hp = Haplotype(chrom=chrom, start=pos1, end=pos2+1, id=hp_id, beta=0.5)
            hp.variants = tuple(
                Variant(start=vr["pos"], end=vr["pos"]+1, id=vr['id'], allele=al)
                for vr, al in zip(snp_vars, allele_combo)
            )
            trans_hp = hp.transform(gts).sum(axis=1)
            combo_ld = np.abs(np.array([pearson_corr_ld(trans_hp, sum_gts[:, snp_id]) for snp_id in combo_id]))
            # we arbitrarily just choose the first value
            bin_idx = np.argmax(combo_ld[0] < intervals)
            if ld_bins[bin_idx] < reps and np.abs(combo_ld[1] - combo_ld[0]) < step:
                log.info(f"Outputting {hp_id} with LD {tuple(combo_ld)} for bin {intervals[bin_idx]}")
                yield hp
                ld_bins[bin_idx] += 1
                count -= 1
                break
        if count == 0:
            break


@click.command()
@click.argument("gts", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-r",
    "--reps",
    type=int,
    default=1,
    show_default=True,
    help="The number of replicates to perform within each LD bin",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="DEBUG",
    show_default=True,
    help="The level of verbosity desired",
)
def main(gts: Path, reps: int = 1):
    """
    Create a column in the GT matrix from the provided haplotype pseudogenotypes

    Output the column to stdout

    \f
    Parameters
    ----------
    gts: Path
        The path to the pgen file containing bialleleic SNP alleles
    """
    log = getLogger("choose", "DEBUG")

    gts = GenotypesPLINK(fname=gts, log=log)
    gts.read()
    gts.check_missing()
    gts.check_biallelic()
    gts.check_phase()

    af = gts.check_maf()
    af_thresh = (af > 0.25) & (af < 0.75)
    gts.subset(variants=gts.variants['id'][af_thresh], inplace=True)

    hps = Haplotypes("/dev/stdout", log=log, haplotype=Haplotype)
    hps.data = {hp.id: hp for hp in find_haps(gts, log, min_ld=0, max_ld=1, reps=reps)}
    hps.write()


if __name__ == "__main__":
    main()
