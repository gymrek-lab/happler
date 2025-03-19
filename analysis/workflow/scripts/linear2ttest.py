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

from happler.tree.assoc_test import NodeResults
from happler.tree.terminator import TTestTerminator
from happler.tree.assoc_test import AssocResults, AssocTestSimple


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
    "-c",
    "--covariance",
    is_flag=True,
    show_default=True,
    default=False,
    help="Use the covariance correction",
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
    hap_id: str = None,
    covariance: bool = False,
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

    log.info("Adjusting betas and stderrs")
    df["beta"] = df["beta"] * linear_stds
    df["stderr"] = df["stderr"] * linear_stds
    parent_df["beta"] = parent_df["beta"] * parent_stds
    parent_df["stderr"] = parent_df["stderr"] * parent_stds

    if covariance:
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
    parent_res = NodeResults.from_np(parent_df.loc[["beta", "pval", "stderr"]])
    results = AssocResults(
        np.array(
            list(df[["beta", "pval", "stderr"]].itertuples(index=False)),
            dtype=AssocTestSimple().return_dtype,
        )
    )
    terminator = TTestTerminator(corrector=None)

    log.info("Computing t-test p-values")
    vals = [
        terminator.compute_val(
            parent_res,
            NodeResults.from_np(val[["beta", "pval", "stderr"]]),
            results,
            idx,
            num_samps,
            num_tests,
            parent_corr=parent_corr,
            short_circuit=False,
        ) for idx, val in df.iterrows()
    ]
    df["pval"] = np.array([(val[0] if val != True else 1) for val in vals])

    if df["pval"].isna().any():
        raise ValueError("Some pvals were NA")

    log.info("Outputting t-test p-values")
    # first, reset column names
    df.rename(columns={v: k for k, v in PLINK_COLS.items()}, inplace=True)
    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()
