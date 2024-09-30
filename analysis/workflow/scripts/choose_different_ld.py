#!/usr/bin/env python
from pathlib import Path
from logging import Logger
from itertools import product

import click
import numpy as np

from haptools.logging import getLogger
from haptools.ld import pearson_corr_ld
from haptools.sim_phenotype import Haplotype
from haptools.data import Genotypes, GenotypesPLINK, Haplotypes, Variant


def find_haps(
    gts: Genotypes,
    log: Logger,
    min_ld: float,
    max_ld: float,
    num_haps: int = 1,
    step: float = 0.1,
    seed: np.random.Generator = np.random.default_rng(None),
):
    log.info(f"Using min_ld {min_ld}, max_ld {max_ld}, num_haps {num_haps}, and step {step}")
    sum_gts = gts.data.sum(axis=2)

    chrom = np.unique(gts.variants["chrom"])
    assert len(chrom) == 1
    chrom = chrom[0]

    intervals = np.arange(max(0, min_ld), min(1, max_ld), step) + step
    ld_bins = {idx: [] for idx in range(len(intervals))}
    count = len(intervals) * num_haps

    combos = np.array(list(range(len(gts.variants))))
    seed.shuffle(combos)
    for combo_id1 in combos:
        combo_id = [combo_id1, combo_id1]
        while combo_id[0] == combo_id[1]:
            combo_id[1] = seed.choice(combos)
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
            if len(ld_bins[bin_idx]) < num_haps:
                log.info(f"Adding {hp_id} with LD {combo_ld} to bin {intervals[bin_idx]}")
                ld_bins[bin_idx].append((combo_ld, hp))
                if len(ld_bins[bin_idx]) == num_haps:
                    yield ld_bins[bin_idx]
                count -= 1
                break
        if count == 0:
            break


@click.command()
@click.argument("gts", type=click.Path(exists=True, path_type=Path))
@click.argument("output", type=click.Path(path_type=Path))
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
    "--num-haps",
    type=int,
    default=1,
    show_default=True,
    help="The number of haplotypes to output",
)
@click.option(
    "--seed",
    type=int,
    default=None,
    show_default="random",
    help="A seed to ensure replicable randomness",
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
    output: Path,
    min_ld: float = 0,
    max_ld: float = 1,
    step: float = 0.1,
    min_af: float = 0.25,
    max_af: float = 0.75,
    num_haps: int = 1,
    seed: int = None,
    verbosity: str = "DEBUG",
):
    """
    Create haplotypes with a range of LD values between their alleles

    LD values are injected into the output from brace expressions
    For example, if we are given a path like "ld_{ld}/happler.hap", then
    we will replace {ld} with the provided ld value
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

    seed = np.random.default_rng(seed)

    for haps in find_haps(
        gts,
        log,
        min_ld=min_ld,
        max_ld=max_ld,
        num_haps=num_haps,
        step=step,
        seed=seed,
    ):
        avg_combo_ld = np.mean([combo_ld for combo_ld, hp in haps])
        out_path = Path(str(output).format(ld=round(avg_combo_ld, 2)))
        out_path.parent.mkdir(parents=True, exist_ok=True)
        hps = Haplotypes(
            out_path, log=log, haplotype=Haplotype,
        )
        hps.data = {hp.id: hp for combo_ld, hp in haps}
        hps.write()


if __name__ == "__main__":
    main()
