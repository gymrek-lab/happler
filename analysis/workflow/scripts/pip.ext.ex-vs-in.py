#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
from haptools.logging import getLogger

import numpy as np
import matplotlib.pyplot as plt

# Read the TSV file into a numpy structured array
data = np.genfromtxt('pips.ext.ex-vs-in.tsv', dtype=None, names=True, delimiter='\t', encoding='utf-8')

# TODO: remake these plots but for different sample size and effect size values

# Sort by delta_bic_train
delta_bic_thresh = 25
data.sort(order='delta_bic_train')

# Calculate percent of points above 0.5 for the bottom row
hap_pip_train_pct_above = np.array([np.sum(data['hap_pip_train'][:x] > 0.5) / len(data['hap_pip_train']) * 100 for x in range(len(data['hap_pip_train']))])
hap_pip_test_pct_above = np.array([np.sum(data['hap_pip_test'][:x] > 0.5) / len(data['hap_pip_test']) * 100 for x in range(len(data['hap_pip_test']))])
# hap_pip_train_pct_above = data['hap_pip_train']
# hap_pip_test_pct_above = data['hap_pip_test']

# Get unique sample sizes and betas
unique_sampsizes = np.unique(data['sampsize'])
unique_betas = np.unique(data['beta'])

# Create a colormap for each sample size
colormaps = {
    sampsize: plt.cm.get_cmap(cmap_name)
    for sampsize, cmap_name in zip(
        unique_sampsizes,
        ['Blues', 'Oranges', 'Reds', 'Purples', 'Greens'][:len(unique_sampsizes)]
    )
}

# Create a function to plot points grouped by color (sample size and beta)
def plot_by_groups(ax, x_data, y_data, sampsizes, betas):
    for sampsize in unique_sampsizes:
        for beta in unique_betas:
            # Find indices for this sampsize/beta combination
            mask = (sampsizes == sampsize) & (betas == beta)
            if np.any(mask):
                # Normalize beta to 0-1 range for colormap
                beta_norm = (beta - np.min(unique_betas)) / (np.max(unique_betas) - np.min(unique_betas))
                cmap = colormaps[sampsize]
                color = cmap(0.3 + 0.7 * beta_norm)

                # Plot all points with this color at once
                ax.scatter(x_data[mask], y_data[mask], color=color, alpha=0.6, s=50)

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
plt.savefig('pips.ext.ex-vs-in.pdf', dpi=300, bbox_inches='tight')
