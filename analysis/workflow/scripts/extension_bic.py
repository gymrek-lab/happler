#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
import numpy as np
import pandas as pd
from haptools import data
from haptools.logging import getLogger
from haptools.ld import pearson_corr_ld

from happler.tree import TreeBuilder
from happler.tree.variant import Variant
from happler.tree.assoc_test import NodeResultsExtra
from happler.tree.terminator import BICTerminator, TTestTerminator
from happler.tree.haplotypes import Haplotype, HapplerHaplotype, HapplerVariant
from happler.tree.assoc_test import AssocResults, AssocTestSimpleSM, AssocTestSimpleSMTScore


def var_gts(
    variant: Variant,
    allele: int,
    gts: data.Genotypes,
):
    return gts.subset(variants=(variant.id,)).data[:,0] == allele

def get_extension_bf_parent(
    hp: data.Haplotype,
    hap_gts: data.Genotypes,
    og_gts: data.Genotypes,
    phen: data.Phenotypes,
    mode: str = "parent-bic",
    maf: float = None,
    log: Logger = None,
):
    if mode != "parent-bic":
        raise ValueError("Unsupported mode")

    log.info("Setting up tree builder")
    hap_tree = TreeBuilder(
        og_gts,
        phen,
        maf=maf,
        method=AssocTestSimpleSM(with_bic=True),
        terminator=BICTerminator(bf_thresh=-float("inf"), log=log),
        indep_thresh=-float("inf"),
        ld_prune_thresh=0.95,
        covariance_correction=False,
        log=log,
    )

    log.info("Locating best delta BIC value among four haplotypes")
    best_bic_value = {"ext_allele": None, "results": None, "bf_val": -float("inf"), "new_haps": None}
    parent_hap = Haplotype(num_samples=len(hap_tree.gens.samples))
    parent_res = None

    # find the variant-allele pairs that give the best haplotype
    vals = hap_tree._find_split_rigid(parent_hap, parent_res)

    if vals is not None:

        for variant, allele, results in vals:
            if variant is None:
                # there were no significant variants!
                continue
            # create a new Haplotype with the variant-allele pair added
            variant_gts = hap_tree.gens.data[:, variant.idx, :2] == allele
            new_parent_hap = parent_hap.append(variant, allele, variant_gts)
            new_vals = hap_tree._find_split_rigid(new_parent_hap, results)
            parent_res = results

            if new_vals is None:
                continue

            for new_variant, new_allele, new_results in new_vals:
                if new_variant is None:
                    continue

                new_variant_gts = hap_tree.gens.data[:, new_variant.idx, :2] == new_allele
                new_hap = new_parent_hap.append(new_variant, new_allele, new_variant_gts)

                num_tests = 1
                parent_corr = 0
                num_samps = int(len(og_gts.samples))
                node_results = new_results
                ext_allele = (new_variant, new_allele, new_results)
                assoc_results = AssocResults(
                    np.array(
                        [(node_results.beta,node_results.pval,node_results.stderr,node_results.bic),],
                        dtype=[
                            ("beta", np.float64),
                            ("pval", object),
                            ("stderr", np.float64),
                            ("bic", np.float64),
                        ],
                    )
                )

                # now, get the BF
                val = BICTerminator().compute_val(
                    parent_res,
                    node_results,
                    assoc_results,
                    0,
                    num_samps,
                    num_tests,
                    parent_corr=parent_corr,
                    short_circuit=False,
                )
                if val != True:
                    bf_val = val[mode.endswith("bic")]
                else:
                    bf_val = float("inf") if mode.endswith("bic") else 0

                if np.isnan(bf_val):
                    raise ValueError("Some BFs were NA")

                if bf_val > best_bic_value["bf_val"]:

                    new_haps = data.Haplotypes(
                        fname=None, haplotype=HapplerHaplotype, variant=HapplerVariant, log=log
                    )
                    new_haps.data = {}
                    hap_node_results = (parent_res.bic, node_results.bic)
                    new_haps.data[hp.id] = HapplerHaplotype.from_happler_haplotype(
                        new_hap, og_gts, hp.id, hap_node_results,
                    )
                    new_haps.data[hp.id].beta = node_results.beta
                    new_haps.data[hp.id].pval = -np.log10(node_results.pval)

                    best_bic_value = {
                        "ext_allele": ext_allele,
                        "results": assoc_results,
                        "bf_val": bf_val,
                        "new_haps": new_haps,
                    }

    return tuple(best_bic_value.values())

def get_extension_bf(
    hp: data.Haplotype,
    hap_gts: data.Genotypes,
    og_gts: data.Genotypes,
    phen: data.Phenotypes,
    mode: str = "bic",
    maf: float = None,
    log: Logger = None,
):
    log.info("Setting up tree builder")
    hap_tree = TreeBuilder(
        og_gts,
        phen,
        maf=maf,
        terminator=BICTerminator(bf_thresh=-float("inf"), log=log),
        indep_thresh=-float("inf"),
        ld_prune_thresh=0.95,
        covariance_correction=False,
        log=log,
    )

    log.info("Setting up delta BIC test")
    # parent node model: y ~ h_hap
    parent = Haplotype.from_haptools_haplotype(hp, og_gts)
    hap_gts_data = hap_gts.data
    if mode.endswith("bic"):
        parent_res = NodeResultsExtra.from_np(
            AssocTestSimpleSM(with_bic=True).run(
                hap_gts_data.sum(axis=2),
                phen.data[:, 0],
            ).data[0]
        )
    elif mode == "tscore":
        parent_res = NodeResultsExtra.from_np(
            AssocTestSimpleSMTScore(with_bic=True).run(
                hap_gts_data.sum(axis=2),
                phen.data[:, 0],
            ).data[0]
        )
    else:
        raise ValueError("Unsupported mode")

    log.info("Running tree builder for a single node")
    ext_allele = list(filter(
        lambda x: x[0] is not None,
        hap_tree._find_split_rigid(parent, parent_res)
    ))
    if len(ext_allele) > 1:
        # if both alleles were unterminated, we choose the one with the best BIC
        ext_allele = min(ext_allele, key=lambda x: x[2].bic)
    else:
        ext_allele = ext_allele[0]
    
    log.info("Obtaining BIC for best allele extension")
    new_allele_gts = og_gts.subset(variants=(ext_allele[0].id,)).data[:, 0]
    if ext_allele[1] == 0:
        new_allele_gts = ~new_allele_gts
    new_hap = parent.append(ext_allele[0], ext_allele[1], new_allele_gts)
    # current node model: y ~ h_hap' where hap' is hap extended by the next best allele
    if mode.endswith("bic"):
        terminator = BICTerminator()
        results = AssocTestSimpleSM(with_bic=True).run(
            new_hap.data.sum(axis=1)[:, np.newaxis],
            phen.data[:, 0],
        )
    elif mode == "tscore":
        terminator = TTestTerminator()
        results = AssocTestSimpleSMTScore(with_bic=True).run(
            new_hap.data.sum(axis=1)[:, np.newaxis],
            phen.data[:, 0],
        )
    else:
        raise ValueError("Unsupported mode")

    log.info("Computing BF values")
    num_tests = 1
    parent_corr = 0
    num_samps = int(len(og_gts.samples))
    node_res = NodeResultsExtra
    node_results = node_res.from_np(results.data[0])
    # check that we were able to recapitulate the results object properly
    if isinstance(ext_allele, NodeResultsExtra):
        assert node_results.beta == ext_allele[2].beta
        assert node_results.stderr == ext_allele[2].stderr
    np.testing.assert_almost_equal(node_results.bic, ext_allele[2].bic, decimal=3)

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
        bf_val = val[mode.endswith("bic")]
    else:
        bf_val = float("inf") if mode.endswith("bic") else 0

    if np.isnan(bf_val):
        raise ValueError("Some BFs were NA")

    new_haps = data.Haplotypes(
        fname=None, haplotype=HapplerHaplotype, variant=HapplerVariant, log=log
    )
    new_haps.data = {}
    # NOTE: this code is broken for --mode bic and needs to be fixed! an extra bic
    # value must be prepended to the following tuple
    hap_node_results = (parent_res.bic, node_results.bic)
    new_haps.data[hp.id] = HapplerHaplotype.from_happler_haplotype(
        new_hap, og_gts, hp.id, hap_node_results,
    )
    new_haps.data[hp.id].beta = node_results.beta
    new_haps.data[hp.id].pval = -np.log10(node_results.pval)

    return ext_allele, results, bf_val, new_haps


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
    "-m",
    "--mode",
    type=click.Choice(["tscore", "bic", "parent-bic"]),
    default="bic",
    show_default=True,
    help="The type of values to compute",
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
    mode: str = "bic",
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Determine the difference in BIC between the hap and its best extension
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
    hap_gts.check_maf(threshold=maf)
    assert len(hap_gts.variants) == 1
    assert phen.samples == og_gts.samples and phen.samples == hap_gts.samples

    # call method to compute BIC
    if mode == "parent-bic":
        ext_allele, results, bf_val, new_haps = get_extension_bf_parent(
            hp, hap_gts, og_gts, phen, mode, maf, log,
        )
    else:
        ext_allele, results, bf_val, new_haps = get_extension_bf(
            hp, hap_gts, og_gts, phen, mode, maf, log,
        )

    log.info("Outputting new .hap file")
    new_haps.fname = output.with_suffix(".hap")
    new_haps.write()

    log.info("Outputting BF values")
    PLINK_COLS = {
        "#CHROM": hp.chrom,
        "POS": ext_allele[0].pos,
        "ID": ext_allele[0].id,
        "OBS_CT": int(len(og_gts.samples)),
        "BETA": results.data["beta"][0],
        "SE": results.data["stderr"][0],
        "P": bf_val,
    }
    df = pd.DataFrame([PLINK_COLS])
    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()
