#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
import numpy as np
import pandas as pd
from haptools import data
from haptools.logging import getLogger

from happler.tree import TreeBuilder
from happler.tree.variant import Variant
from happler.tree.haplotypes import Haplotype
from happler.tree.assoc_test import NodeResultsExtra
from happler.tree.terminator import BICTerminator
from happler.tree.assoc_test import AssocTestSimpleSM

@click.command()
@click.argument("hap", type=click.Path(exists=True, path_type=Path))
@click.argument("hap_gts", type=click.Path(exists=True, path_type=Path))
@click.argument("og_gts", type=click.Path(exists=True, path_type=Path))
@click.argument("phenotype", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--maf",
    type=float,
    default=None,
    show_default="no filtering",
    help="Ignore variants with a MAF below this threshold",
)
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
    hap: Path,
    hap_gts: Path,
    og_gts: Path,
    phenotype: Path,
    maf: float = None,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Determine the difference in BIC (as a BF) between the hap and its best extension
    """
    log = getLogger("extension-bic", verbosity)

    log.info("Loading haplotypes, phenotypes, and genotypes")
    hp = list(data.Haplotypes.load(hap).data.values())[0]
    phen = data.Phenotypes.load(phenotype)
    og_gts = data.GenotypesPLINK(og_gts)
    og_gts.read(samples=set(phen.samples))
    phen.subset(samples=og_gts.samples, inplace=True)
    og_gts.check_missing(discard_also=True)
    og_gts.check_biallelic(discard_also=True)
    og_gts.check_maf(threshold=maf, discard_also=True)
    og_gts.check_phase()
    hap_gts = data.GenotypesPLINK.load(hap_gts)
    hap_gts.check_missing()
    hap_gts.check_biallelic()
    hap_gts.check_maf(threshold=maf)
    hap_gts.check_phase()
    assert len(hap_gts.variants) == 1
    assert phen.samples == og_gts.samples and phen.samples == hap_gts.samples

    log.info("Setting up delta BIC test")
    # parent node model: y ~ h_hap
    parent = Haplotype.from_haptools_haplotype(hp, og_gts)
    parent_res = NodeResultsExtra.from_np(
        AssocTestSimpleSM(with_bic=True).run(
            hap_gts.data.sum(axis=2),
            phen.data[:, 0],
        ).data[0]
    )

    log.info("Setting up tree builder")
    hap_tree = TreeBuilder(
        og_gts,
        phen,
        maf=maf,
        terminator=BICTerminator(bf_thresh=0, log=log),
        indep_thresh=1,
        ld_prune_thresh=0.95,
        covariance_correction=False,
        log=log,
    )

    log.info("Running tree builder for a single node")
    ext_allele = list(filter(
        lambda x: x[0] is not None,
        hap_tree._find_split_rigid(parent, parent_res)
    ))
    if len(ext_allele) > 1:
        # if both alleles were unterminated, we choose the one with the best BIC
        ext_allele = max(ext_allele, key=lambda x: x[2].bic)
    else:
        ext_allele = ext_allele[0]
    
    log.info("Obtaining BIC for best allele extension")
    new_allele_gts = og_gts.subset(variants=(ext_allele[0].id,)).data[:, 0]
    if ext_allele[1] == 0:
        new_allele_gts = ~new_allele_gts
    new_hap = parent.append(ext_allele[0], ext_allele[1], new_allele_gts)
    # current node model: y ~ h_hap' where hap' is hap extended by the next best allele
    results = AssocTestSimpleSM(with_bic=True).run(
        new_hap.data.sum(axis=1)[:, np.newaxis],
        phen.data[:, 0],
    )

    log.info("Computing BF values")
    num_tests = 1
    parent_corr = 0
    num_samps = int(len(og_gts.samples))
    node_res = NodeResultsExtra
    node_results = node_res.from_np(results.data[0])
    terminator = BICTerminator()
    # check that we were able to recapitulate the results object properly
    assert node_results == ext_allele[2]

    # now, get the BF
    val = terminator.compute_val(
        parent_res,
        node_results,
        results,
        0,
        num_samps,
        num_tests,
        parent_corr=parent_corr,
        short_circuit=False,
    )
    if val != True:
        bf_val = val[1]
    else:
        bf_val = float("inf")

    if np.isnan(bf_val):
        raise ValueError("Some BFs were NA")

    log.info("Outputting t-test p-values")
    PLINK_COLS = {
        "#CHROM": hp.chrom,
        "POS": ext_allele[0].pos,
        "ID": ext_allele[0].id,
        "OBS_CT": num_samps,
        "BETA": results.data["beta"][0],
        "SE": results.data["stderr"][0],
        "P": bf_val,
    }
    df = pd.DataFrame([PLINK_COLS])
    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()
