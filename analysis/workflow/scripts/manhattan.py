#!/usr/bin/env python

import sys
import click
import numpy as np
import pandas as pd
from pathlib import Path
import statsmodels.api as sm
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype


@click.command()
@click.option(
    "-l",
    "--linear",
    multiple=True,
    type=click.Path(exists=True, path_type=Path),
    default=[Path("-")],
    show_default="stdin",
    help="PLINK2 .linear files containing the assocation results",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("-"),
    show_default="stdout",
    help="A PNG file containing a Manhattan plot of the results",
)
@click.option(
    "-i",
    "--id",
    "ids",
    multiple=True,
    type=str,
    default=tuple(),
    show_default="no IDs",
    help="Which variant IDs should we highlight in red?",
)
@click.option(
    "-b",
    "--orange-id",
    "orange_ids",
    multiple=True,
    type=str,
    default=tuple(),
    show_default="no IDs",
    help="Which variant IDs should we highlight in blue?",
)
@click.option(
    "--label/--no-label",
    is_flag=True,
    default=True,
    show_default=True,
    help="Whether to label the points by their IDs, as well",
)
def main(linear=[sys.stdin], output=sys.stdout, ids=tuple(), orange_ids=tuple(), label=True):
    """
    Create a manhattan plot from the results of a PLINK2 GWAS
    """
    plink_cols = {
        "#CHROM": "chromosome",
        "POS": "pos",
        "ID": "id",
        "REF": "ref",
        "ALT": "alt",
        "A1": "allele",
        "TEST": "test",
        "OBS_CT": "num_samps",
        "BETA": "beta",
        "SE": "se",
        "T_STAT": "tstat",
        "P": "pval",
        "ERRCODE": "error",
    }
    keep_cols = ["chromosome", "pos", "id", "beta", "se", "pval"]

    # create the plot
    fig, ax = plt.subplots(1, len(linear), figsize=(14, 8), sharex=True, sharey=True)

    if len(linear) == 1:
        ax = (ax,)

    # parse the .linear files
    dfs = {}
    max_pval = -1
    red_ids = ids
    for idx, linear_fname in enumerate(linear):
        df = pd.read_csv(
            linear_fname,
            sep="\t",
            header=0,
            names=plink_cols.values(),
            usecols=keep_cols,
        ).sort_values("pos")
        pos_range = max(df["pos"]) - min(df["pos"])
        label_distance = pos_range/45
        # replace NaN with inf
        df["pval"] = df["pval"].fillna(np.inf)
        df['-log10(p)'] = -np.log10(df["pval"])
        # replace -infinity values with 0
        df['-log10(p)'].replace([-np.inf], 0, inplace=True)
        df.chromosome = df.chromosome.astype('category')
        df.chromosome = df.chromosome.astype(
            CategoricalDtype(sorted(map(int, df.chromosome.dtype.categories)), ordered=True)
        )
        df = df.sort_values('chromosome')
        # create the plot using pandas and add it to the figure
        df[~df["id"].isin(red_ids + orange_ids)].plot(kind='scatter', x='pos', y='-log10(p)', ax=ax[idx])
        # plot red ids if there are any
        if red_ids:
            v_ids = df[df["id"].isin(red_ids)]['id']
            x_ids = df[df["id"].isin(red_ids)]['pos']
            y_ids = df[df["id"].isin(red_ids)]['-log10(p)']
            if np.any(np.isinf(y_ids)):
                raise ValueError(f"The p-values for {red_ids} are too powerful!")
            ax[idx].scatter(x_ids, y_ids, color='red', marker='o', s=20)
            if label:
                for v_id, x_id, y_id in zip(v_ids, x_ids, y_ids):
                    ax[idx].annotate(v_id, (x_id+label_distance, y_id))
        # plot blue ids if there are any
        if orange_ids:
            v_ids = df[df["id"].isin(orange_ids)]['id']
            x_ids = df[df["id"].isin(orange_ids)]['pos']
            y_ids = df[df["id"].isin(orange_ids)]['-log10(p)']
            if np.any(np.isinf(y_ids)):
                raise ValueError(f"The p-values for {orange_ids} are too powerful!")
            ax[idx].scatter(x_ids, y_ids, color='orange', marker='o', s=20)
            if label:
                for v_id, x_id, y_id in zip(v_ids, x_ids, y_ids):
                    ax[idx].annotate(v_id, (x_id+label_distance, y_id))
        # set title and perform cleanup/secondary tasks
        ax_name = Path(Path(linear_fname.stem).stem).stem
        ax[idx].title.set_text(ax_name.replace('-', " + "))
        ax[idx].set(xlabel=None, ylabel=None)
        dfs[ax_name] = df
        max_val = ax[idx].get_ylim()[1]
        if max_pval < max_val:
            max_pval = max_val
    df = pd.concat(dfs, ignore_index=True)

    # set the y-axis limit so that both axes have the same limit
    for idx, linear_fname in enumerate(linear):
        ax[idx].set_ylim(top=max_pval)

    # save the graph
    fig.supxlabel('Chromosomal Position')
    fig.supylabel('-log10(p value)')
    plt.savefig(output)


if __name__ == "__main__":
    main()
