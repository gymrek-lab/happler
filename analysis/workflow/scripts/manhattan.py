#!/usr/bin/env python

import sys
import click
import numpy as np
import pandas as pd
from pathlib import Path
import statsmodels.api as sm
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype


AXIS_FONTSIZE = 6
TITLE_FONTSIZE = 5
AXIS_LABELSIZE = 2.5
LABEL_FONTSIZE = 4
TICK_FONTSIZE = 4
POINT_SIZE = 0.75


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
    "--orange-Hids",
    is_flag=True,
    default=False,
    show_default=True,
    help="Add any IDs of the form 'H:digit:' to the list of orange IDs",
)
@click.option(
    "--label/--no-label",
    is_flag=True,
    default=True,
    show_default=True,
    help="Whether to label the points by their IDs, as well",
)
@click.option(
    "--small",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to shrink the plot to a smaller size",
)
@click.option(
    "-t",
    "--titles",
    multiple=True,
    type=str,
    default=tuple(),
    show_default="infer from IDs",
    help="Which titles should be given to each plot?",
)
@click.option(
    "-a",
    "--alpha",
    type=float,
    default=None,
    show_default="do not draw line",
    help="alpha threshold to depict on graph",
)
def main(
    linear=[sys.stdin],
    output=sys.stdout,
    ids=tuple(),
    orange_ids=tuple(),
    orange_hids=False,
    label=True,
    small=False,
    titles=tuple(),
    alpha=None,
):
    """
    Create a manhattan plot from the results of a PLINK2 GWAS
    """
    plink_cols = {
        "#CHROM": "chromosome",
        "POS": "pos",
        "ID": "id",
        "P": "pval",
    }
    keep_cols = ["#CHROM", "POS", "ID", "P"]
    # create the plot
    fig, ax = plt.subplots(
        1, len(linear), sharex=True, sharey=True, constrained_layout=True,
    )

    # parse the .linear files
    dfs = {}
    max_pval = -1
    red_ids = ids
    for idx, linear_fname in enumerate(linear):
    
        # first, handle cases where there may be more than one .linear file
        cur_ax = ax
        if cur_ax is tuple:
            cur_ax = cur_ax[idx]

        df = pd.read_csv(
            linear_fname,
            sep="\t",
            header=0,
            usecols=keep_cols,
        ).rename(columns=plink_cols)
        df = df.sort_values("pos")
        pos_range = max(df["pos"]) - min(df["pos"])
        label_distance = pos_range/17
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
        # retrieve any haplotype IDs (of the form 'H:digit:')
        if orange_hids:
            orange_ids += tuple(df.id[df.id.str.contains("^H\d+$")])
        # create the plot using pandas and add it to the figure
        if small:
            df[~df["id"].isin(red_ids + orange_ids)].plot(
                kind='scatter', x='pos', y='-log10(p)', ax=cur_ax, s=POINT_SIZE,
            )
        else:
            df[~df["id"].isin(red_ids + orange_ids)].plot(
                kind='scatter', x='pos', y='-log10(p)', ax=cur_ax,
            )
        # plot red ids if there are any
        if red_ids:
            v_ids = df[df["id"].isin(red_ids)]['id']
            x_ids = df[df["id"].isin(red_ids)]['pos']
            y_ids = df[df["id"].isin(red_ids)]['-log10(p)']
            if np.any(np.isinf(y_ids)):
                raise ValueError(f"The p-values for {red_ids} are too powerful!")
            if small:
                cur_ax.scatter(x_ids, y_ids, color='red', marker='o', s=POINT_SIZE)
            else:
                cur_ax.scatter(x_ids, y_ids, color='red', marker='o', s=20)
            if label:
                for v_id, x_id, y_id in zip(v_ids, x_ids, y_ids):
                    if small:
                        cur_ax.annotate(v_id, (x_id+label_distance, y_id), fontsize=LABEL_FONTSIZE)
                    else:
                        cur_ax.annotate(v_id, (x_id+label_distance, y_id))
        # plot blue ids if there are any
        if orange_ids:
            v_ids = df[df["id"].isin(orange_ids)]['id']
            x_ids = df[df["id"].isin(orange_ids)]['pos']
            y_ids = df[df["id"].isin(orange_ids)]['-log10(p)']
            if np.any(np.isinf(y_ids)):
                raise ValueError(f"The p-values for {orange_ids} are too powerful!")
            if small:
                cur_ax.scatter(x_ids, y_ids, color='orange', marker='o', s=POINT_SIZE)
            else:
                cur_ax.scatter(x_ids, y_ids, color='orange', marker='o', s=20)
            if label:
                for v_id, x_id, y_id in zip(v_ids, x_ids, y_ids):
                    if small:
                        cur_ax.annotate(v_id, (x_id+label_distance, y_id), fontsize=LABEL_FONTSIZE)
                    else:
                        cur_ax.annotate(v_id, (x_id+label_distance, y_id))
        # set title and perform cleanup/secondary tasks
        ax_name = Path(Path(linear_fname.stem).stem).stem
        if ax_name == "haplotype":
            ax_name = "Haplotype effect"
        if titles:
            print(f"using title {titles[idx]}")
            ax_name = titles[idx]
        if small:
            cur_ax.set_title(ax_name.replace('-', " + "), fontdict={
                'fontsize': TITLE_FONTSIZE
            })
            cur_ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
            cur_ax.tick_params(axis='both', which='minor', labelsize=TICK_FONTSIZE)
            cur_ax.xaxis.get_offset_text().set_fontsize(TICK_FONTSIZE)
        else:
            cur_ax.set_title(ax_name.replace('-', " + "))
        cur_ax.set(xlabel=None, ylabel=None)
        dfs[ax_name] = df
        max_val = cur_ax.get_ylim()[1]
        if max_pval < max_val:
            max_pval = max_val
    df = pd.concat(dfs, ignore_index=True)

    # if we haven't already flipped alpha to the log scale, do it now
    if alpha is not None and alpha > 0 and alpha < 0.5:
        alpha = -np.log10(alpha)

    # set the y-axis limit so that both axes have the same limit
    for idx, linear_fname in enumerate(linear):
        cur_ax.set_ylim(top=max_pval)
        if alpha is not None:
            cur_ax.axhline(y=alpha, color='red')

    # save the graph
    if small:
        fig.supxlabel('Chromosomal Position', fontsize=AXIS_FONTSIZE)
        fig.supylabel('$-log_{10} P-value$', fontsize=AXIS_FONTSIZE)
    else:
        fig.supxlabel('Chromosomal Position')
        fig.supylabel('$-log_{10} P-value$')
    if small:
        fig.set_size_inches(2.65, 2.2)
        plt.savefig(output, bbox_inches="tight", pad_inches=0.03)
    else:
        fig.set_size_inches(5, 3.33)
        plt.savefig(output)


if __name__ == "__main__":
    main()
