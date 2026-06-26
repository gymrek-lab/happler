#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
from haptools.logging import getLogger

import numpy as np
import matplotlib.pyplot as plt


@click.command()
@click.argument("pips_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-s",
    "--sampsize",
    type=int,
    default=None,
    show_default=True,
    help="Subset to points with this sample size",
)
@click.option(
    "-b",
    "--beta",
    type=float,
    default=None,
    show_default=True,
    help="Subset to points with this beta",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,
    show_default="input.pdf",
    help="A PNG or PDF file containing the desired plots",
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
    pips_file: Path,
    sampsize: int = None,
    beta: int = None,
    output: Path = None,
    verbosity: str = "INFO"
):
    """
    Create a heatmap plot that visualizes a haplotype from a .hap file
    """
    log = getLogger("pip.ext", level=verbosity)

    # Read the TSV file into a numpy structured array
    data = np.genfromtxt(str(pips_file), dtype=None, names=True, delimiter='\t', encoding='utf-8')
    if output is None:
        output = pips_file.with_suffix(".pdf")

    # filter to points with the selected sample size and beta
    if sampsize is not None:
        passing = data["sampsize"] == sampsize
        log.info(f"Subsetting to {len(passing)} points with sample size: {sampsize}")
        data = data[passing]
    if beta is not None:
        passing = data["beta"] == beta
        log.info(f"Subsetting to {len(passing)} points with beta: {beta}")
        data = data[passing]

    # Sort by delta_bic_train
    delta_bic_thresh = 25
    data.sort(order='delta_bic_train')

    # Calculate percent of points above 0.5 for the bottom row
    hap_pip_train_pct_above = np.array([np.sum(data['hap_pip_train'][x:] > 0.5) / len(data['hap_pip_train']) * 100 for x in range(len(data['hap_pip_train']))])
    hap_pip_test_pct_above = np.array([np.sum(data['hap_pip_test'][x:] > 0.5) / len(data['hap_pip_test']) * 100 for x in range(len(data['hap_pip_test']))])
    # hap_pip_train_pct_above = data['hap_pip_train']
    # hap_pip_test_pct_above = data['hap_pip_test']

    # Create figure with 2x2 subplots
    fig, axes = plt.subplots(3, 2, figsize=(9.5, 12))

    # Plot 1: include_pip_train vs exclude_pip_train
    axes[0, 0].scatter(data['exclude_pip_train'], data['include_pip_train'], alpha=0.6)
    # plot_by_groups(
    #     axes[0, 0],
    #     data['exclude_pip_train'],
    #     data['include_pip_train'],
    #     data['sampsize'],
    #     data['beta'],
    #     'exclude_pip_train',
    #     'include_pip_train',
    # )
    axes[0, 0].set_xlabel("exclude_pip_train")
    axes[0, 0].set_ylabel("include_pip_train")
    axes[0, 0].grid(True, alpha=0.3)

    # Plot 2: include_pip_test vs exclude_pip_test
    axes[0, 1].scatter(data['exclude_pip_test'], data['include_pip_test'], alpha=0.6)
    # plot_by_groups(
    #     axes[0, 1],
    #     data['exclude_pip_test'],
    #     data['include_pip_test'],
    #     data['sampsize'],
    #     data['beta'],
    #     'exclude_pip_test',
    #     'include_pip_train',
    # )
    axes[0, 1].set_xlabel("exclude_pip_test")
    axes[0, 1].set_ylabel("include_pip_test")
    axes[0, 1].grid(True, alpha=0.3)

    # Plot 3: hap_pip_train vs delta_bic_train
    axes[1, 0].scatter(data['delta_bic_train'], data['hap_pip_train'], alpha=0.6)
    axes[1, 0].set_xlabel('delta_bic_train')
    axes[1, 0].set_ylabel('hap_pip_train_percent_above_0.5')
    axes[1, 0].axvline(x=delta_bic_thresh, color="red")
    axes[1, 0].grid(True, alpha=0.3)

    # Plot 4: hap_pip_test vs delta_bic_train
    axes[1, 1].scatter(data['delta_bic_train'], data['hap_pip_test'], alpha=0.6)
    axes[1, 1].set_xlabel('delta_bic_train')
    axes[1, 1].set_ylabel('hap_pip_test_percent_above_0.5')
    axes[1, 1].axvline(x=delta_bic_thresh, color="red")
    axes[1, 1].grid(True, alpha=0.3)

    # Plot 5: hap_pip_train_percent vs delta_bic_train
    axes[2, 0].plot(data['delta_bic_train'], hap_pip_train_pct_above)
    axes[2, 0].set_xlabel('delta_bic_train')
    axes[2, 0].set_ylabel('hap_pip_train_percent_above_0.5')
    axes[2, 0].axvline(x=delta_bic_thresh, color="red")
    axes[2, 0].grid(True, alpha=0.3)

    # Plot 6: hap_pip_test_percent vs delta_bic_train
    axes[2, 1].plot(data['delta_bic_train'], hap_pip_test_pct_above)
    axes[2, 1].set_xlabel('delta_bic_train')
    axes[2, 1].set_ylabel('hap_pip_test_percent_above_0.5')
    axes[2, 1].axvline(x=delta_bic_thresh, color="red", label=f"threshold: ({delta_bic_thresh})")
    axes[2, 1].grid(True, alpha=0.3)

    # Adjust layout to prevent overlap
    plt.tight_layout()
    fig.legend(loc="lower right")

    # Save and show the plot
    plt.savefig(output, dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    main()
