#!/usr/bin/env python
from typing import Tuple
from pathlib import Path
from logging import Logger
from decimal import Decimal

import click
import numpy as np
import pandas as pd
from haptools import data
from haptools.logging import getLogger
from haptools.ld import pearson_corr_ld

from happler.tree.assoc_test import NodeResults, NodeResultsExtra
from happler.tree.terminator import TTestTerminator, BICTerminator
from happler.tree.assoc_test import (
    AssocResults,
    AssocTestSimple,
    AssocTestSimpleSM,
    AssocTestSimpleCovariates
)


PLINK_COLS = {
    "#CHROM": "chromosome",
    "POS": "pos",
    "ID": "id",
    "OBS_CT": "samples",
    "BETA": "beta",
    "SE": "stderr",
    "P": "pval",
}

def load_linear_file(linear_fname: Path):
    keep_cols = list(PLINK_COLS.keys())
    df = pd.read_csv(
        linear_fname,
        sep="\t",
        header=0,
        usecols=keep_cols,
        converters={"P": lambda val: Decimal(val if val != "NA" else 1)},
    ).rename(columns=PLINK_COLS)
    df = df.sort_values("pos")
    df["pval"] = df["pval"].fillna(np.inf)
    return df


@click.command()
@click.argument("linear", type=click.Path(exists=True, path_type=Path))
@click.argument("parent", type=click.Path(exists=True, path_type=Path))
@click.argument("linear_gts", type=click.Path(exists=True, path_type=Path))
@click.argument("parent_gts", type=click.Path(exists=True, path_type=Path))
@click.argument("phenotype", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-i",
    "--id",
    "hap_id",
    type=str,
    show_default="first",
    help=(
        "The ID of the variant to choose as the parent haplotype"
    ),
)
@click.option(
    "-m",
    "--mode",
    type=click.Choice(["tscore", "covariance", "bic", "interact-bic"]),
    default="tscore",
    show_default=True,
    help="The type of values to compute",
)
@click.option(
    "--child-gts",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    show_default=True,
    help=(
        "The genotypes file for the child haplotype. "
        "Only needed if --mode is interact-bic"
    ),
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
    linear: Path,
    parent: Path,
    linear_gts: Path,
    parent_gts: Path,
    phenotype: Path,
    hap_id: str = None,
    mode: str = "tscore",
    child_gts: Path = None,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Convert p-values from plink2 --glm into t-test p-values using the given parent

    The parent file should come from running plink2 --glm on haptools transform
    """
    log = getLogger("linear2ttest", verbosity)

    log.info("Loading linear files")
    df = load_linear_file(linear)
    parent_df = load_linear_file(parent)
    if hap_id is not None:
        parent_df = parent_df[parent_df.id == hap_id]
    parent_df = parent_df.iloc[0]

    log.info("Loading genotypes")
    parent_gts = data.GenotypesPLINK.load(parent_gts)
    # reorder to match
    parent_gts.subset(variants=(parent_df["id"],), inplace=True)
    parent_stds = parent_gts.data.sum(axis=2).std(axis=0)
    linear_gts = data.GenotypesPLINK.load(linear_gts)
    # reorder to match
    linear_gts.subset(variants=tuple(df["id"]), inplace=True)
    linear_stds = linear_gts.data.sum(axis=2).std(axis=0)
    # also load child gts if needed
    if mode == "interact-bic":
        if child_gts is None:
            raise ValueError("Child genotypes file is required for interact-bic mode")
        child_gts = data.GenotypesPLINK.load(child_gts)
        # reorder to match
        child_gts.subset(variants=tuple(child_gts.variants["id"]), inplace=True)
        child_stds = child_gts.data.sum(axis=2).std(axis=0)

    log.info("Adjusting betas and stderrs")
    df["beta"] = df["beta"] * linear_stds
    df["stderr"] = df["stderr"] * linear_stds
    parent_df["beta"] = parent_df["beta"] * parent_stds
    parent_df["stderr"] = parent_df["stderr"] * parent_stds

    if mode == "covariance":
        log.info("Computing correlations")
        parent_corr = pearson_corr_ld(
            linear_gts.data.sum(axis=2),
            parent_gts.data.sum(axis=2),
        )[:,0]
    else:
        parent_corr = 0

    log.info("Setting up t-tests")
    num_tests = 1
    num_samps = int(parent_df.samples)
    if mode == "bic" or mode == "interact-bic":
        phen = data.Phenotypes.load(phenotype)
        if mode == "bic":
            # parent node model: y ~ h_parent
            parent_res = NodeResultsExtra.from_np(
                    AssocTestSimpleSM(with_bic=True).run(
                    parent_gts.data.sum(axis=2),
                    phen.data[:, 0],
                ).data[0]
            )
            # current node model: y ~ h_hap
            results = AssocTestSimpleSM(with_bic=True).run(
                linear_gts.data.sum(axis=2),
                phen.data[:, 0],
            )
        elif mode == "interact-bic":
            child_covar = child_gts.data.sum(axis=2)
            parent_and_child_covar = data.GenotypesPLINK.merge_variants(
                (parent_gts, child_gts), fname=None
            ).data.sum(axis=2)
            if True:
                # parent node model: y ~ h_parent + z_child
                parent_res = NodeResultsExtra.from_np(
                        AssocTestSimpleCovariates(covars=child_covar, with_bic=True).run(
                        parent_gts.data.sum(axis=2),
                        phen.data[:, 0],
                    ).data[0]
                )
                # current node model: y ~ h_hap
                results = AssocTestSimpleSM(with_bic=True).run(
                    linear_gts.data.sum(axis=2),
                    phen.data[:, 0],
                )
            elif False:
                # parent node model: y ~ h_parent + z_child
                parent_res = NodeResultsExtra.from_np(
                        AssocTestSimpleCovariates(covars=child_covar, with_bic=True).run(
                        parent_gts.data.sum(axis=2),
                        phen.data[:, 0],
                    ).data[0]
                )
                # current node model: y ~ h_hap + h_parent + z_child
                results = AssocTestSimpleCovariates(covars=parent_and_child_covar, with_bic=True).run(
                    linear_gts.data.sum(axis=2),
                    phen.data[:, 0],
                ).data[0]
            elif False:
                # parent node model: y ~ h_hap + h_parent + z_child
                parent_res = NodeResultsExtra.from_np(
                        AssocTestSimpleCovariates(covars=parent_and_child_covar, with_bic=True).run(
                        linear_gts.data.sum(axis=2),
                        phen.data[:, 0],
                    ).data[0]
                )
                # current node model: y ~ h_hap
                results = AssocTestSimpleSM(with_bic=True).run(
                    linear_gts.data.sum(axis=2),
                    phen.data[:, 0],
                )
        node_res = NodeResultsExtra
        terminator = BICTerminator()
    else:
        parent_res = NodeResults.from_np(parent_df.loc[["beta", "pval", "stderr"]])
        results = AssocResults(
            np.array(
                list(df[["beta", "pval", "stderr"]].itertuples(index=False)),
                dtype=AssocTestSimple().return_dtype,
            )
        )
        node_res = NodeResults
        terminator = TTestTerminator(corrector=None)

    log.info("Computing t-test p-values")
    vals = [
        terminator.compute_val(
            parent_res,
            node_res.from_np(val),
            results,
            idx,
            num_samps,
            num_tests,
            parent_corr=parent_corr,
            short_circuit=False,
        ) for idx, val in enumerate(results.data)
    ]
    bic = int(mode == "bic" or mode == "interact-bic")
    df["pval"] = np.array([(val[bic] if val != True else 1) for val in vals])

    if df["pval"].isna().any():
        raise ValueError("Some pvals were NA")

    log.info("Outputting t-test p-values")
    # first, reset column names
    df.rename(columns={v: k for k, v in PLINK_COLS.items()}, inplace=True)
    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()
