#!/usr/bin/env python
import logging
import warnings
from pathlib import Path
from decimal import Decimal

import click
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from haptools.logging import getLogger

from snakemake_io import glob_wildcards, remove_regexes


def get_metrics(
    metrics_file: Path,
    log: logging.Logger = None
) -> float:
    """
    Extract the metrics from the metrics files

    Parameters
    ----------
    metrics_file: Path
        The path to a .hap file containing a set of haplotypes
    log: logging.Logger, optional
        A logging object to write any debugging and error messages

    Returns
    -------
    npt.NDArray
        Mixed dtype array with named columns
    """
    with open(metrics_file) as mf:
        column_names = mf.readlines()[0].strip().split("\t")
    return np.loadtxt(
        metrics_file,
        skiprows=1,
        delimiter="\t",
        dtype=[(n, np.float64) for n in column_names],
    )


@click.command()
@click.argument("metrics_files", type=click.Path(path_type=Path))
@click.option(
    "--use-flex-axes-limits",
    is_flag=True,
    show_default=True,
    default=False,
    help=(
        "Allow for flexible axes limits. Makes it harder to compare across plots but"
        "ensures that all points are visible"
    ),
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A .png file containing the output plot",
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
    metrics_files: Path,
    use_flex_axes_limits: bool = False,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Plot metrics.tsv files from midway_manhattan_summary.py in simple plots

    Each file is inferred from brace expressions injected into the paths to the files
    For example, a path like "{region}.{type}/{tswitch}/out.{rep}.{name}.glm.linear"
    will infer the wildcards and match them between cases. Wildcards will be expanded
    across all inputs. So, for example, if the case_type is "tswitch", then these two
    files will be matched together:
    {region}.{type}/1/out.{rep}.{name}.tsv
    {region}.{type}/2/out.{rep}.{name}.tsv
    """
    log = getLogger("midway_summary_metrics", level=verbosity)

    # extract parameters and parameter values by globbing wildcards
    # params will be a dictionary mapping parameter names to lists of values
    params = dict(glob_wildcards(metrics_files)._asdict())
    dtypes = {k: "U30" for k in params.keys()}
    # ensure ints and floats are intepreted appropriately so that sorting works later
    if "sampsize" in dtypes:
        dtypes["sampsize"] = np.int64
    if "beta" in dtypes:
        dtypes["beta"] = np.float16
    log.info(f"Extracted paramter values {tuple(dtypes.keys())}")
    # convert the dictionary to a numpy mixed dtype array
    params = np.array(list(zip(*params.values())), dtype=list(dtypes.items()))
    params.sort()
    # convert back to U30 so that file matching works
    params = params.astype([(k, "U30") for k in dtypes.keys()])

    metrics_files_wo_regexes = remove_regexes(str(metrics_files))
    get_fname = lambda path, param_set: Path(str(path).format(**dict(zip(dtypes.keys(), param_set))))


    log.info(f"Extracting metrics from {len(params)} metrics files")
    metrics = np.array([
        get_metrics(
            get_fname(metrics_files_wo_regexes, params[idx]),
            log=log
        )
        for idx in range(len(params))
    ])

    log.info("Creating plots")
    fig, axs = plt.subplots(
        nrows=1, ncols=5,
        sharex=True, figsize=((7/3)*5, 3),
        constrained_layout=True, tight_layout=False,
    )

    log.info("Writing to axes")
    betas = np.unique(params["beta"])
    for beta in betas:
        params_beta = params[params["beta"] == beta]
        metrics_beta = metrics[params["beta"] == beta]

        axs[0].plot(params_beta["sampsize"], metrics_beta["FPR"], "o", label=beta)
        if not use_flex_axes_limits:
            axs[0].set_ylim((-0.001, 0.051))
        axs[0].set_ylabel("False positive rate")

        axs[1].plot(params_beta["sampsize"], metrics_beta["Recall"], "o")
        if not use_flex_axes_limits:
            axs[1].set_ylim((-0.001, 1.001))
        axs[1].set_ylabel("Recall")

        axs[2].plot(params_beta["sampsize"], metrics_beta["Significance Threshold"], "o")
        max_alpha = max(metrics_beta["Significance Threshold"])
        if not use_flex_axes_limits:
            if max_alpha > 1:
                axs[2].set_ylim((-2.005, 10.005))
            elif max_alpha < 0.25:
                axs[2].set_ylim((-0.005, 0.255))
            else:
                axs[2].set_ylim((-0.005, 1.005))
        axs[2].set_ylabel("Significance threshold")

        axs[3].plot(params_beta["sampsize"], metrics_beta["AUROC"], "o")
        if not use_flex_axes_limits:
            axs[3].set_ylim((0.895, 1.005))
        axs[3].set_ylabel("AUROC")

        axs[4].plot(params_beta["sampsize"], metrics_beta["Average Precision"], "o")
        if not use_flex_axes_limits:
            axs[4].set_ylim((0.895, 1.005))
        axs[4].set_ylabel("Average Precision")

    log.info("Writing figure")
    fig.supxlabel("Sample size")
    fig.legend(loc="lower right", ncol=len(betas), fontsize="small")
    plt.savefig(output)


if __name__ == "__main__":
    main()
