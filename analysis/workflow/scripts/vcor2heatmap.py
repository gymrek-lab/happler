#!/usr/bin/env python
from pathlib import Path

import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@click.command()
@click.argument("vcor_file", type=click.Path(exists=True, path_type=Path))
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
    vcor_file: Path,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Make an LD correlation table plot of a vcor file

    Example usage
    -------------
    log_file="outs/out-geuvadis-newthresholds/11_333589-1333589/happler/run/ENSG00000177697.12/logs/run"
    variants="$(grep -F '[label="' "$log_file" | grep -v root | sed 's+^. \[label="++;s+\\n{.*$++' | sort -u)"
    plink2 --extract <(echo "$variants") --ld-window-kb 9999999 --ld-window-r2 0 --out pairwise_hap.ld --pfile data/geuvadis/geuvadis_ensemble_phasing.maf0.1 --r2-unphased
    workflow/scripts/vcor2heatmap.py -o pairwise_hap.ld.png pairwise_hap.ld.vcor
    """
    vcor_df = pd.read_csv(vcor_file, sep="\t")
    variants = pd.concat([vcor_df['ID_A'], vcor_df['ID_B']]).unique()
    variant_map = {variant: i for i, variant in enumerate(variants)}
    n = len(variants)
    corr_matrix = np.eye(n)
    for _, row in vcor_df.iterrows():
        i = variant_map[row['ID_A']]
        j = variant_map[row['ID_B']]
        corr = row['UNPHASED_R2']
        corr_matrix[i, j] = corr
        corr_matrix[j, i] = corr
    fig, ax = plt.subplots(figsize=(10, 10))
    heatmap = ax.matshow(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1)
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(variants, rotation=90, fontsize=10)
    ax.set_yticklabels(variants, fontsize=10)
    ax.tick_params(axis='x', bottom=True, top=False, labelbottom=True, labeltop=False)
    ax.set_title('Variant Correlation Heatmap', pad=20)
    cbar = plt.colorbar(heatmap)
    cbar.set_label('Correlation (r)', rotation=270, labelpad=15)
    plt.tight_layout()
    plt.savefig(output)


if __name__ == "__main__":
    main()
