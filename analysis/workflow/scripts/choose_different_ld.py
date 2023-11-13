#!/usr/bin/env python
import random
from pathlib import Path
from logging import Logger
from itertools import combinations, product

import click
import numpy as np

from haptools.logging import getLogger
from haptools.ld import pearson_corr_ld
from haptools.sim_phenotype import Haplotype
from haptools.data import Genotypes, GenotypesPLINK, Haplotypes, Variant


def random_combinations(iterable, r):
    """
        Random selection from itertools.combinations(iterable, r)
        Copied from the ex recipe at https://docs.python.org/library/itertools.html
    """
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)


def find_haps(
    gts: Genotypes,
    log: Logger,
    min_ld: float,
    max_ld: float,
    reps: int = 1,
    step: float = 0.1,
):
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
    for combo_id1 in combos:
        combo_id = [combo_id1, combo_id1]
        while combo_id[0] == combo_id[1]:
            combo_id[1] = np.random.choice(combos)
        snp_vars = gts.variants[combo_id]
        hp_id = "-".join(snp_vars['id'])
        log.debug(f"Trying {hp_id}")
        pos1, pos2 = snp_vars["pos"]
        for allele_combo in product(*snp_vars["alleles"]):
            hp = Haplotype(chrom=chrom, start=pos1, end=pos2+1, id=hp_id, beta=0.99)
            hp.variants = tuple(
                Variant(start=vr["pos"], end=vr["pos"]+1, id=vr['id'], allele=al)
                for vr, al in zip(snp_vars, allele_combo)
            )
            combo_ld = np.abs(pearson_corr_ld(*tuple(sum_gts[:, snp_id] for snp_id in combo_id)))
            # which bin does this LD value fall in?
            bin_idx = np.argmax(combo_ld < intervals)
            if ld_bins[bin_idx] < reps:
                log.info(f"Outputting {hp_id} with LD {combo_ld} for bin {intervals[bin_idx]}")
                yield combo_ld, hp
                ld_bins[bin_idx] += 1
                count -= 1
                break
        if count == 0:
            break


@click.command()
@click.argument("gts", type=click.Path(exists=True, path_type=Path))
@click.argument("output_dir", type=click.Path(path_type=Path))
@click.option(
    "-r",
    "--reps",
    type=int,
    default=1,
    show_default=True,
    help="The number of replicates to perform within each LD bin",
)
@click.option(
    "--min-ld",
    type=float,
    default=0,
    show_default=True,
    help="The minimum LD value to allow",
)
@click.option(
    "--max-ld",
    type=float,
    default=1,
    show_default=True,
    help="The maximum LD value to allow",
)
@click.option(
    "--step",
    type=float,
    default=0.1,
    show_default=True,
    help="The step size between each LD bin",
)
@click.option(
    "--min-af",
    type=float,
    default=0.25,
    show_default=True,
    help="The minimum LD value to allow",
)
@click.option(
    "--max-af",
    type=float,
    default=0.75,
    show_default=True,
    help="The maximum LD value to allow",
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
    gts: Path,
    output_dir: Path,
    reps: int = 1,
    min_ld: float = 0,
    max_ld: float = 1,
    step: float = 0.1,
    min_af: float = 0.25,
    max_af: float = 0.75,
    verbosity: str = "DEBUG",
):
    """
    Create haplotypes with a range of LD values between their alleles
    """
    log = getLogger("choose", verbosity)

    gts = GenotypesPLINK(fname=gts, log=log)
    gts.read()
    gts.check_missing()
    gts.check_biallelic()
    gts.check_phase()

    af = gts.check_maf()
    af_thresh = (af > min_af) & (af < max_af)
    gts.subset(variants=gts.variants['id'][af_thresh], inplace=True)

    # create output directory
    output_dir.mkdir()

    for combo_ld, hp in find_haps(
        gts,
        log,
        min_ld=min_ld,
        max_ld=max_ld,
        reps=reps,
        step=step,
    ):
        out_dir = output_dir / f"ld_{round(combo_ld, 2)}"
        out_dir.mkdir()
        hps = Haplotypes(
            out_dir / "haplotype.hap", log=log, haplotype=Haplotype,
        )
        hps.data = {hp.id: hp}
        hps.write()


if __name__ == "__main__":
    main()
