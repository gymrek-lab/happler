#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
from haptools.logging import getLogger

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import ConnectionPatch


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
    "--include-raw-scatter",
    is_flag=True,
    show_default=True,
    default=False,
    help="Should we include a scatterplot containing the raw hap PIPs?",
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
    include_raw_scatter: bool = False,
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
    delta_bic_thresh = 20
    data.sort(order='delta_bic_train')

    # Calculate percent of points above 0.5 for the bottom row
    hap_pip_train_pct_above = np.array([np.sum(data['hap_pip_train'][x:] > 0.5) / len(data['hap_pip_train']) * 100 for x in range(len(data['hap_pip_train']))])
    hap_pip_test_pct_above = np.array([np.sum(data['hap_pip_test'][x:] > 0.5) / len(data['hap_pip_test']) * 100 for x in range(len(data['hap_pip_test']))])
    # hap_pip_train_pct_above = data['hap_pip_train']
    # hap_pip_test_pct_above = data['hap_pip_test']

    # Create figure with 2x2 subplots
    # If include_raw_scatter is true, then we add one more row with two axes
    fig, axes = plt.subplots(
        2+int(include_raw_scatter), # height
        2, # width
        figsize=(
            4.75*2, # width
            4*(2+int(include_raw_scatter)), # height
        ),
    )

    if axes.ndim <= 1:
        axes = axes[np.newaxis, :]

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
    axes[0, 0].set_xlabel("PIP of causal parent SNP when haplotype is excluded")
    axes[0, 0].set_ylabel("PIP of causal parent SNP when haplotype is included")
    axes[0, 0].set_title("Training Cohort")
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
    axes[0, 1].set_xlabel("PIP of causal parent SNP when haplotype is excluded")
    axes[0, 1].set_ylabel("PIP of causal parent SNP when haplotype is included")
    axes[0, 1].set_title("Validation Cohort")
    axes[0, 1].grid(True, alpha=0.3)

    if include_raw_scatter:
        # Plot 3: hap_pip_train vs delta_bic_train
        axes[1, 0].scatter(data['delta_bic_train'], data['hap_pip_train'], alpha=0.6)
        axes[1, 0].set_xlabel('Delta BIC in training cohort')
        axes[1, 0].set_ylabel('Non-causal haplotype PIPs')
        axes[1, 0].axvline(x=delta_bic_thresh, color="red")
        axes[1, 0].grid(True, alpha=0.3)

        # Plot 4: hap_pip_test vs delta_bic_train
        axes[1, 1].scatter(data['delta_bic_train'], data['hap_pip_test'], alpha=0.6)
        axes[1, 1].set_xlabel('Delta BIC in training cohort')
        axes[1, 1].set_ylabel('Non-causal haplotype PIPs')
        axes[1, 1].axvline(x=delta_bic_thresh, color="red")
        axes[1, 1].grid(True, alpha=0.3)

    # Plot 5: hap_pip_train_percent vs delta_bic_train (zoomed out)
    axes[1+int(include_raw_scatter), 0].plot(data['delta_bic_train'], hap_pip_test_pct_above, color="orange")
    axes[1+int(include_raw_scatter), 0].plot(data['delta_bic_train'], hap_pip_train_pct_above, color="blue")
    axes[1+int(include_raw_scatter), 0].set_xlabel('Delta BIC in training cohort')
    axes[1+int(include_raw_scatter), 0].set_ylabel('Percent of non-causal haplotype PIPs > 0.5')
    axes[1+int(include_raw_scatter), 0].axvline(x=delta_bic_thresh, color="red", linestyle="--")
    axes[1+int(include_raw_scatter), 0].grid(True, alpha=0.3)

    # Plot 6: hap_pip_test_percent vs delta_bic_train (zoomed in)
    axes[1+int(include_raw_scatter), 1].plot(data['delta_bic_train'], hap_pip_test_pct_above, label="Validation cohort", color="orange")
    axes[1+int(include_raw_scatter), 1].set_xlabel('Delta BIC in training cohort')
    #axes[1+int(include_raw_scatter), 1].set_ylabel('Percent of non-causal haplotype PIPs > 0.5')
    axes[1+int(include_raw_scatter), 1].grid(True, alpha=0.3)
    axes[1+int(include_raw_scatter), 1].set_ylim(bottom=-0.05)
    minor_tick_val = axes[1+int(include_raw_scatter), 1].get_yticks().max()
    axes[1+int(include_raw_scatter), 1].plot(data['delta_bic_train'], hap_pip_train_pct_above, label="Training cohort", color="blue")
    axes[1+int(include_raw_scatter), 1].set_ylim(top=minor_tick_val)
    axes[1+int(include_raw_scatter), 1].axvline(x=delta_bic_thresh, color="red", label=f"Delta BIC threshold: {delta_bic_thresh}", linestyle="--")
    axes[1+int(include_raw_scatter), 1].legend(loc="upper right")

    # add y-axis tick mark to the hap_pip_train plot at the max tick mark of the hap_pip_test plot
    # axes[1+int(include_raw_scatter), 0].yaxis.set_minor_locator(ticker.FixedLocator([minor_tick_val]))

    ax1 = axes[1+int(include_raw_scatter), 0]
    ax2 = axes[1+int(include_raw_scatter), 1]
    con1 = ConnectionPatch(
        xyA=(0, minor_tick_val), coordsA=ax1.get_yaxis_transform(),
        xyB=(0, 1.0), coordsB=ax2.transAxes,
        axesA=ax1, axesB=ax2,
        color="gray",
        linestyle="--",
        linewidth=1.5,
    )
    con2 = ConnectionPatch(
        xyA=(0, -0.05), coordsA=ax1.get_yaxis_transform(),
        xyB=(0, 0), coordsB=ax2.transAxes,
        axesA=ax1, axesB=ax2,
        color="gray",
        linestyle="--",
        linewidth=1.5,
    )
    fig.add_artist(con1)
    fig.add_artist(con2)

    # add panel labels
    panel_labels = ["A", "B", "C", "D", "E", "F"]
    visible_axes = [ax for row in axes for ax in row if ax.has_data()]
    for label, ax in zip(panel_labels, visible_axes[:4]):
        ax.text(
            -0.14, 1.05, label,
            transform=ax.transAxes,
            fontsize=20,
            fontweight="bold",
            va="top",
            ha="right",
        )

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save and show the plot
    plt.savefig(output, dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    main()
