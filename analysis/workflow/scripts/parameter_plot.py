#!/usr/bin/env python
import copy
import pickle
import warnings
from typing import Union
from pathlib import Path
from logging import Logger
from dataclasses import replace

import click
import matplotlib
import numpy as np
matplotlib.use('Agg')
import numpy.typing as npt
from scipy.stats import sem
from haptools import logging
import matplotlib.pyplot as plt
from haptools.ld import pearson_corr_ld
from numpy.lib import recfunctions as rfn
from scipy.optimize import linear_sum_assignment
from haptools.data import Genotypes, GenotypesVCF, GenotypesPLINK, Haplotypes, Phenotypes

from happler.tree.haplotypes import (
    HapplerVariant as Variant,
    HapplerHaplotype as Haplotype,
)

from snakemake_io import glob_wildcards
from extension_bic import get_extension_bf


DTYPES = {
    "beta": np.float64,
    "ld": np.float64,
    "rep": np.uint8,
    "gt": "U5",
    "alpha": np.uint32,
    "samp": np.uint32,
    "num_haps": np.uint8,
}

# LOG_SCALE = {"alpha",}
LOG_SCALE = {}

plt.rcParams['figure.dpi'] = 400  # Set the figure DPI to 300
plt.rcParams['savefig.dpi'] = plt.rcParams['figure.dpi']  # Set the DPI for saving figures
warnings.filterwarnings("ignore", message="Mean of empty slice", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="All-NaN axis encountered", category=RuntimeWarning)

def match_haps(
    gts: Genotypes,
    observed: Union[Haplotype, Haplotypes],
    causal: Union[Haplotype, Haplotypes],
) -> tuple:
    """
    Match the observed and causal haplotypes to maximize best pairwise LD

    This function uses scipy.optimize.linear_sum_assignment to find the optimal
    match between causal and observed haplotypes to maximize the sum of the pairwise LD

    Parameters
    ----------
    gts: Genotypes
        The set of genotypes for each of the variants in the haplotypes
    observed: Haplotype|Haplotypes
        A path to a .hap file containing haplotypes output by happler
    causal: Haplotype|Haplotypes
        A path to a .hap file containing a simulated causal haplotype

    Returns
    -------
    npt.NDArray
        An array of the LD values for each pair of observed and causal haplotypes
    npt.NDArray
        An array of the row indices indicating the optimal match of causal haplotypes
        to observed haplotypes
    npt.NDArray
        An array of the column indices indicating the optimal match of causal haplotypes
        to observed haplotypes
    npt.NDArray
        The index of the causal haplotype that each observed haplotype is matched to
    npt.NDArray
        An array of boolean values indicating whether each observed haplotype has a match
        in the causal haplotypes
    """
    if isinstance(observed, Haplotypes):
        obs = observed.transform(gts).data.sum(axis=2)
    else:
        obs = observed.transform(gts).sum(axis=1)[:,np.newaxis]
    if isinstance(causal, Haplotypes):
        exp = causal.transform(gts).data.sum(axis=2)
    else:
        exp = causal.transform(gts).sum(axis=1)[:,np.newaxis]
    # Compute the pairwise LD between each observed and causal haplotype
    # The rows and columns of ld_mat both correspond to the observed and causal haps in
    # the order they were given. For example, if there is 1 observed hap and 3 causal
    # haps, then there should be four rows and four columns in ld_mat
    # Retrieve only the correlation of observed haps (as rows) vs causal haps (as cols)
    ld_mat = np.abs(pearson_corr_ld(obs, exp))
    if exp.shape[1] == 1:
        row_idx, col_idx = np.argmax(ld_mat, axis=0), np.array([0], dtype=np.int64)
        assert (row_idx, col_idx) == linear_sum_assignment(ld_mat, maximize=True)
    else:
        # Return the ld_mat idxs of the optimal match of each observed hap to a causal hap
        row_idx, col_idx = linear_sum_assignment(ld_mat, maximize=True)
    # now deal with the extras (observed haps with no causal hap match)
    extras_labels = np.zeros(ld_mat.shape[0], dtype=np.uint8)
    extras_labels_bool = np.zeros(ld_mat.shape[0], dtype=bool)
    missing_rows = np.setdiff1d(range(ld_mat.shape[0]), row_idx)
    extras_labels[row_idx] = col_idx
    extras_labels_bool[row_idx] = True
    # If there are more observed haps than causal haps, then we just assign the
    # remaining haps to the best causal hap that matches for each of them
    row_idx = np.concatenate((row_idx, missing_rows))
    missing_cols = ld_mat[missing_rows].argmax(axis=1)
    col_idx = np.concatenate((col_idx, missing_cols))
    extras_labels[missing_rows] = missing_cols
    extras_labels_bool[missing_rows] = False
    return ld_mat, row_idx, col_idx, extras_labels, extras_labels_bool


def hap_subsets(hap: Haplotype, log: Logger = None):
    """
    Yield subsets of a haplotype that increase in size by one variant/allele at a time

    Parameters
    ----------
    hap: Haplotype
    log: Logger, optional

    Yields
    ------
    Iterator[Haplotype, float]
        The haplotype and its Bayes factor value compared to its parent
    """
    variants = hap.variants
    for i in range(1, len(variants)):
        hp = replace(hap)
        hp.variants = variants[:i]
        bf = variants[i-1].score - variants[i].score
        yield hp, bf


def optimize_all_observed_haps(
    gts: Genotypes,
    observed: Haplotypes,
    causal: Haplotype,
    log: Logger = None,
):
    optimized_haps = []
    best_bfs = []

    # Pre-transform causal hap once
    causal_vec = causal.transform(gts).sum(axis=1)[:, np.newaxis]

    for hap_id, hap in observed.data.items():
        best_ld = 0
        best_sub = hap
        best_bf = None

        for hp_sub, bf_sub in hap_subsets(hap):
            sub_obs = hp_sub.transform(gts).sum(axis=1)[:, np.newaxis]
            ld = np.abs(pearson_corr_ld(sub_obs, causal_vec))[0, 0]

            if ld > best_ld:
                best_ld = ld
                best_sub = hp_sub
                best_bf = bf_sub

        optimized_haps.append(best_sub)
        best_bfs.append(best_bf if best_bf is not None else 0.0)

    return optimized_haps, np.array(best_bfs)


def get_best_ld(
    gts: GenotypesVCF,
    observed_hap: Path,
    causal_hap: Path,
    region: str = None,
    maf: float = None,
    observed_id: str = None,
    causal_id: str = None,
    pheno: Path = None,
    log: Logger = None
):
    """
    Compute the best LD between the observed haplotypes and the causal haplotype

    Parameters
    ----------
    gts: GenotypesVCF
        A GenotypesVCF (or GenotypesPLINK) object from which to load genotypes
    observed_hap: Path
        A path to a .hap file containing haplotypes output by happler
    causal_hap: Path
        A path to a .hap file containing a simulated causal haplotype
    region: str, optional
        Only load genotypes from this region (or the entire file, otherwise)
    maf: float, optional
        Ignore variants with an MAF below this threshold
    observed_id: str, optional
        The ID to load from the observed_hap file (or just all of the IDs, otherwise)
    causal_id: str, optional
        The ID to load from the causal_hap file (or just the first hap, otherwise)
    pheno: Path, optional
        Should we also output the optimal bayes factor thresholds for the observed and
        causal haps? This assumes that the V scores in the hap files are BIC values and
        only works if the causal haps files have a single hap in them.
    log: Logger, optional
        A logging module to pass to haptools

    Returns
    -------
    observed_ld : np.ndarray
        The best LD value for each observed haplotype vs the causal
    extras_labels : np.ndarray
        The index of the causal haplotype that each observed haplotype is matched to
    labels_bool : np.ndarray
        Boolean array indicating which observed haps have a causal match
    best_thresholds : np.ndarray, optional
        The Bayes factor threshold that gave the best LD (if optimize_thresh is True)
    """
    # load the causal haplotype given by 'causal_id' or just the first hap
    causal_hap = Haplotypes(causal_hap, log=log)
    causal_hap.read(haplotypes=(set((causal_id,)) if causal_id is not None else None))

    # load the observed haplotype given by 'observed_id' or just all of the haps
    observed_hap = Haplotypes(observed_hap, haplotype=Haplotype, variant=Variant, log=log)
    observed_hap.read(haplotypes=(set((observed_id,)) if observed_id is not None else None))

    # if there are no observed haplotypes, then we can't compute LD
    if len(observed_hap.data) == 0:
        return (
            np.array([np.nan,]),
            np.array([np.iinfo(np.uint8).max,], dtype=np.uint8),
            np.array([False,], dtype=bool),
        )

    # load the genotypes
    variants = {
        v.id for hap in causal_hap.data
        for v in causal_hap.data[hap].variants
    }
    variants.update({
        v.id for hap in observed_hap.data
        for v in observed_hap.data[hap].variants
    })
    gts._var_idx = None
    og_gts = copy.deepcopy(gts)
    gts.read(variants=variants, region=region)
    gts.check_phase()
    gts.check_missing()
    gts.check_biallelic()

    # if requested, try to optimize the LD values by incrementally subsetting the
    # haplotypes that matched with causal haplotypes
    if pheno is not None:
        assert len(causal_hap.data) == 1, "Only one causal haplotype is supported when --phenos is set."
        causal_hap_single = next(iter(causal_hap.data.values()))
        hap_gts = causal_hap.transform(gts)
        phen = Phenotypes.load(pheno)
        og_gts.read(samples=set(phen.samples))
        phen.subset(samples=og_gts.samples, inplace=True)
        og_gts.check_missing(discard_also=True)
        og_gts.check_biallelic(discard_also=True)
        og_gts.check_maf(threshold=maf, discard_also=True)
        og_gts.check_phase()
        _, _, causal_best_bf = get_extension_bf(
            causal_hap_single, hap_gts, og_gts, phen, mode="bic", maf=maf, log=log,
        )

        # Optimize all observed haps
        # Take subsets of them until we find the subset that 
        optimized_hap_list, best_bfs = optimize_all_observed_haps(
            gts, observed_hap, causal_hap_single, log=log
        )

        # Replace observed haplotype data with optimized ones
        for hap_id, opt_hap in zip(observed_hap.data.keys(), optimized_hap_list):
            observed_hap.data[hap_id] = opt_hap

    # compute LD between every observed haplotype and the causal haplotype
    observed_ld, best_row_idx, best_col_idx, extras_labels, labels_bool = match_haps(
        gts, observed_hap, causal_hap,
    )

    # get the strongest possible LD for each observed hap with a causal hap
    # note: incidentally, extras_labels[best_row_idx] will always be the same as best_col_idx
    observed_ld = observed_ld[best_row_idx, best_col_idx]
    extra_labels = extras_labels[best_row_idx]
    labels_bool = labels_bool[best_row_idx]

    if pheno is not None:
        # reorder to match observed_ld order
        best_bfs = best_bfs[best_row_idx]
        # mark unmatched haps
        best_bfs[~labels_bool] = np.nan
        return observed_ld, extra_labels, labels_bool, best_bfs, causal_best_bf
    else:
        return observed_ld, extra_labels, labels_bool


def get_finemap_metrics(
    metrics_path: Path,
    keep_hap_ids: bool = False,
    log: Logger = None
):
    """
    Parse data from a TSV containing metrics from fine-mapping

    Parameters
    ----------
    metrics: Path
        The path to a metrics.tsv file
    keep_hap_ids: bool, optional
        If True, do not remove the hap IDs in the first column
    log: Logger, optional
        A logging module to pass to haptools
    """
    # The metrics are:
    # 1) What is the ID of the hap?
    # 2) What is the PIP of the hap?
    # 3) Does the observed/causal hap get the highest PIP?
    # 4) What is the next best PIP in the credible set, excluding the hap?
    # 5) Is the observed/causal hap in a credible set? If so, what is it's index?
    # 6) How many credible sets are there?
    # 7) What is the purity of the credible set with the observed/causal hap?
    # 8) What is the length of the credible set with the observed/causal hap?

    # If you add or remove any metrics here, make sure to also change the
    # "num_expected_vals" in get_metrics_mean_std so it knows the number of metrics
    dtype = [
        ("hap_id", "U30"),
        ("pip", np.float64),
        ("has_highest_pip", np.bool_),
        ("best_variant_pip", np.float64),
        ("in_credible_set", np.bool_),
        ("num_credible_sets", np.uint8),
        ("cs_purity", np.float64),
        ("cs_length", np.float64),
    ]
    null_val = np.array([
        (np.nan, np.nan, False, np.nan, False, 0, np.nan, 0),
    ], dtype=dtype)
    # if there are no observed haplotypes, then we can't compute LD
    if not metrics_path.exists():
        return null_val
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", UserWarning)
        metrics = np.atleast_1d(np.loadtxt(fname=metrics_path, delimiter=" ", dtype=dtype))
    if not metrics.shape[0]:
        return null_val
    if not keep_hap_ids:
        # Remove the 'hap_id' column (index 0)
        metrics = rfn.drop_fields(metrics, 'hap_id')
    return metrics


def get_metrics_mean_std(metrics_vals, num_expected_vals: int = 8):
    """
    Compute the mean and standard error of the metrics

    Parameters
    ----------
    metrics_vals: npt.NDArray
        A numpy array containing the metrics for each observed hap
    num_expected_vals: int, optional
        The number of expected metrics

    Returns
    -------
    npt.NDArray
        A numpy array containing the mean of the metrics
    npt.NDArray
        A numpy array containing the standard error of the metrics
    """
    if len(metrics_vals) == 0:
        n = num_expected_vals
        return np.array([np.nan,]*n), np.array([np.nan,]*n)
    nan_vals = np.isnan(metrics_vals["pip"])
    metrics_vals = metrics_vals[~nan_vals]
    mean_vals = np.array(
        [np.nanmean(metrics_vals[field]) for field in metrics_vals.dtype.names],
    )
    sem_vals = np.array(
        [sem(metrics_vals[field], nan_policy="omit") for field in metrics_vals.dtype.names],
    )
    return mean_vals, sem_vals


def count_shared(observed, causal, log, observed_id: str = None, causal_id: str = None):
    observed = Haplotypes(observed, log=log)
    observed.read(haplotypes=(set((observed_id,)) if observed_id is not None else None))
    obs_vars = set(
        var.id
        for hap in observed.data
        for var in observed.data[hap].variants
    )

    causal = Haplotypes(causal, log=log)
    causal.read(haplotypes=(set((causal_id,)) if causal_id is not None else None))
    causal_vars = set(
        var.id
        for hap in causal.data
        for var in causal.data[hap].variants
    )

    # return the number of shared SNPs and the number of expected SNPs
    return len(obs_vars & causal_vars), len(causal_vars)


def draw_vertical_gridlines(fig, axs, num_ticks):
    for x in range(num_ticks):
        color = "#eeeeee" if x % 2 == 0 else "#cccccc"
        for ax in axs:
            ax.axvline(
                x,
                color=color,
                linewidth=0.5,
                zorder=0,
            )


def remove_same_valued_columns(arr):
    """
    remove any columns that have the same value across all rows
    """
    cols_to_keep = []
    for name in arr.dtype.names:
        col = arr[name]
        if not np.all(col == col[0]):
            cols_to_keep.append(name)

    # Create new array with only non-constant columns
    new_dtype = [(name, arr.dtype[name]) for name in cols_to_keep]
    filtered_arr = np.empty(arr.shape, dtype=new_dtype)
    for name in cols_to_keep:
        filtered_arr[name] = arr[name]

    return filtered_arr


def plot_params(
    params,
    vals,
    vals_sem,
    val_title: str,
    val_color,
    metrics: dict = None,
    hide_extras: bool = False,
    best_bf_vals=None,
):
    """
    Plot vals against parameter values

    Parameters
    ----------
    params: npt.NDArray
        A numpy array of mixed dtype. Each column contains the values of a parameter
    vals: npt.NDArray
        A numpy array of numpy arrays containing the values of the plot
    vals_sem: npt.NDArray
        A numpy array of numpy arrays containing the standard error of the values
    val_title: str
        The name of the value that we are plotting
    val_color: npt.NDArray
        Whether to color each point corresponding to this value
    metrics: dict, optional
        A dict mapping fine-mapping metrics to tuples of two lists of numpy arrays
        containing means and standard errors of fine-mapping metrics for each observed hap
    hide_extras: bool, optional
        Whether to hide the observed haps that have no causal hap match
    best_bf_vals: npt.NDArray
        A numpy array of numpy arrays containing the best BF values
    """
    figsize = matplotlib.rcParams["figure.figsize"]
    params = remove_same_valued_columns(params)
    if len(params) > 75:
        figsize[0] = len(params) / 15
    num_rows = len(params.dtype)
    if metrics is not None:
        num_rows += len(metrics)
    if best_bf_vals is not None:
        num_rows += 1
    if num_rows > 5:
        figsize[1] = num_rows * 1.25

    fig, axs = plt.subplots(
        nrows=num_rows + 1, ncols=1,
        sharex=True, figsize=figsize,
    )
    fig.subplots_adjust(hspace=0)

    # main LD values
    x_vals = [j for j, arr in enumerate(vals) for i in arr]
    val_color = ["green" if j else "red" for i in val_color for j in i]
    vals = np.concatenate(vals)
    vals_sem = np.concatenate(vals_sem)
    for v in sorted(zip(x_vals, vals, vals_sem, val_color), key=lambda x: x[3], reverse=True):
        if hide_extras and v[3] != "green":
            continue
        axs[0].errorbar(v[0], v[1], yerr=v[2], marker="o", c=v[3], markersize=3)
    axs[0].set_ylabel(val_title, rotation="horizontal", ha="right")
    axs[0].set_xticks(range(0, len(np.unique(x_vals))))
    axs[0].set_xticklabels([])
    axs[0].set_ylim(None, 1)

    # parameter plots
    for idx, param in enumerate(params.dtype.names):
        val_title = param
        if param in LOG_SCALE:
            params[param] = -np.log10(params[param])
            val_title = "-log " + val_title
        axs[idx+1].plot(params[param], "-")
        axs[idx+1].set_ylabel(val_title, rotation="horizontal", ha="right")

    offset = len(params.dtype) + 1

    # fine-mapping metrics
    if metrics is not None:
        for idx, metric in enumerate(metrics.keys()):
            curr_ax = axs[offset + idx]
            vals = np.concatenate(metrics[metric][0])
            vals_sem = np.concatenate(metrics[metric][1])
            for v in sorted(zip(x_vals, vals, vals_sem, val_color), key=lambda x: x[3], reverse=True):
                if hide_extras and v[3] != "green":
                    continue
                curr_ax.errorbar(v[0], v[1], yerr=v[2], marker="o", c=v[3], markersize=3)
            curr_ax.set_ylabel(metric, rotation="horizontal", ha="right")

    # best_bf panel
    if best_bf_vals is not None:
        curr_ax = axs[-1]
        vals = np.concatenate(best_bf_vals)
        vals_sem = np.full_like(vals, np.nan)
        for v in sorted(zip(x_vals, vals, vals_sem, val_color), key=lambda x: x[3], reverse=True):
            if hide_extras and v[3] != "green":
                continue
            curr_ax.errorbar(v[0], v[1], yerr=None, marker="o", c=v[3], markersize=3)
        curr_ax.set_ylabel("Best BF", rotation="horizontal", ha="right")

    draw_vertical_gridlines(fig, axs, len(np.unique(x_vals)))
    return fig


def plot_params_simple(
    params,
    vals,
    vals_sem,
    val_title: str,
    val_color,
    hide_extras: bool = False,
):
    """
    Plot vals against only a single column of parameter values

    Parameters
    ----------
    params: npt.NDArray
        A numpy array containing the values of a parameter
    vals: npt.NDArray
        A numpy array of numpy arrays containing the values of the plot
    vals_sem: npt.NDArray
        A numpy array of numpy arrays containing the standard error of the values
    val_title: str
        The name of the value that we are plotting
    val_color: npt.NDArray
        Whether to color each point corresponding to this value
    hide_extras: bool, optional
        Whether to hide the observed haps that have no causal hap match
    """
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(3.5, 3))
    params = remove_same_valued_columns(params)
    val_xlabel = params.dtype.names[0]
    if val_xlabel in LOG_SCALE:
        params = -np.log10(params)
        val_xlabel = "-log " + val_xlabel
    params = [j for j, arr in zip(params, vals) for i in arr]
    val_color = ["green" if j else "red" for i in val_color for j in i]
    vals = np.concatenate(vals)
    vals_sem = np.concatenate(vals_sem)
    # plot params against vals
    for v in sorted(
        zip(params, vals, vals_sem, val_color),
        key=lambda x: x[3],
        reverse=True
    ):
        if hide_extras and v[3] != "green":
            continue
        for v_0 in v[0]:
            ax.errorbar(v_0, v[1], yerr=v[2], marker="o", c=v[3], markersize=5)
    ax.set_ylabel(val_title, color="g")
    ax.set_xlabel(val_xlabel)
    ax.set_ylim(0, 1.03)
    return fig


def group_by_rep(
    params,
    vals,
    causal_idxs,
    bools,
    metrics=None,
    extra_vals=None,
    extra_vals_2=None,
):
    """
    Group replicates with identical parameter values

    Compute mean and standard error of the values for each group

    Parameters
    ----------
    params: npt.NDArray
        A numpy array of mixed dtype. Each column contains the values of a parameter
    vals: npt.NDArray
        A numpy array containing the values of the plot
    causal_idxs: npt.NDArray
        The index of the causal haplotype that each observed haplotype belongs to
    bools: npt.NDArray
        An array of boolean values indicating whether each observed haplotype has a match
    metrics: npt.NDArray
        A numpy array of numpy arrays containing the metrics for each observed hap
    """
    other_param_names = [name for name in params.dtype.names if name != "rep"]
    grouped_params = np.unique(params[other_param_names])
    get_mean_std = lambda x: (np.nanmean(x), sem(x, nan_policy="omit") if len(x) > 1 else np.nan)
    filter_causal_idxs = lambda x: np.unique(x[x != np.iinfo(np.uint8).max])

    grouped_vals = []
    grouped_sem = []
    grouped_bools = []
    grouped_metrics = []
    grouped_metrics_sem = []
    grouped_extra = [] if extra_vals is not None else None
    grouped_extra_2 = [] if extra_vals_2 is not None else None

    def safe_index(val, indices, mask):
        """Ensure val and mask shapes match before indexing."""
        val = np.asarray(val)
        indices = np.asarray(indices)
        mask = np.asarray(mask)
        if val.shape != mask.shape:
            if val.size == 1:
                val = np.full(mask.shape, val.item())
            else:
                return np.array([], dtype=val.dtype)
        return val[mask]

    for group in grouped_params:
        subgrouped_vals = []
        subgrouped_sem = []
        subgrouped_bools = []
        subgrouped_metrics = []
        subgrouped_metrics_sem = []
        subgrouped_extra = []
        subgrouped_extra_2 = []

        curr_vals = vals[params[other_param_names] == group]
        curr_causal_idxs = causal_idxs[params[other_param_names] == group]
        curr_bools = bools[params[other_param_names] == group]
        curr_metrics = metrics[params[other_param_names] == group] if metrics is not None else None
        curr_extra_vals = extra_vals[params[other_param_names] == group] if extra_vals is not None else None
        curr_extra_vals_2 = extra_vals_2[params[other_param_names] == group] if extra_vals_2 is not None else None

        for causal_idx in filter_causal_idxs(np.concatenate(curr_causal_idxs)):
            try:
                curr_vals_causal = np.concatenate([
                    safe_index(val, indices, (indices == causal_idx) & curr_bool)
                    for val, indices, curr_bool in zip(curr_vals, curr_causal_idxs, curr_bools)
                ])
            except ValueError:
                continue
            val_mean, val_sem = get_mean_std(curr_vals_causal)
            subgrouped_vals.append(val_mean)
            subgrouped_sem.append(val_sem)
            subgrouped_bools.append(True)

            if curr_metrics is not None:
                curr_metrics_causal = np.concatenate([
                    safe_index(val, indices, (indices == causal_idx) & curr_bool)
                    for val, indices, curr_bool in zip(curr_metrics, curr_causal_idxs, curr_bools)
                ])
                metrics_mean, metrics_sem = get_metrics_mean_std(curr_metrics_causal)
                subgrouped_metrics.append(metrics_mean)
                subgrouped_metrics_sem.append(metrics_sem)

            if curr_extra_vals is not None:
                curr_extra_causal = np.concatenate([
                    safe_index(val, indices, (indices == causal_idx) & curr_bool)
                    for val, indices, curr_bool in zip(curr_extra_vals, curr_causal_idxs, curr_bools)
                ])
                subgrouped_extra.append(np.nanmean(curr_extra_causal))

            if curr_extra_vals_2 is not None:
                curr_extra_causal_2 = np.concatenate([
                    safe_index(val, indices, (indices == causal_idx) & curr_bool)
                    for val, indices, curr_bool in zip(curr_extra_vals_2, curr_causal_idxs, curr_bools)
                ])
                subgrouped_extra_2.append(np.nanmean(curr_extra_causal_2))

            # unmatched haps
            for unmatched_idx in range(max([(~i).sum() for i in curr_bools])):
                try:
                    curr_vals_unmatched = np.array([
                        safe_index(val, indices, (indices == causal_idx) & ~curr_bool)[unmatched_idx]
                        for val, indices, curr_bool in zip(curr_vals, curr_causal_idxs, curr_bools)
                        if unmatched_idx < safe_index(val, indices, (indices == causal_idx) & ~curr_bool).shape[0]
                    ])
                except ValueError:
                    continue
                val_mean, val_sem = get_mean_std(curr_vals_unmatched)
                subgrouped_vals.append(val_mean)
                subgrouped_sem.append(val_sem)
                subgrouped_bools.append(False)

                if curr_metrics is not None:
                    curr_metrics_unmatched = np.array([
                        safe_index(val, indices, (indices == causal_idx) & ~curr_bool)[unmatched_idx]
                        for val, indices, curr_bool in zip(curr_metrics, curr_causal_idxs, curr_bools)
                        if unmatched_idx < safe_index(val, indices, (indices == causal_idx) & ~curr_bool).shape[0]
                    ])
                    metrics_mean, metrics_sem = get_metrics_mean_std(curr_metrics_unmatched)
                    subgrouped_metrics.append(metrics_mean)
                    subgrouped_metrics_sem.append(metrics_sem)

                if curr_extra_vals is not None:
                    curr_extra_unmatched = np.array([
                        safe_index(val, indices, (indices == causal_idx) & ~curr_bool)[unmatched_idx]
                        for val, indices, curr_bool in zip(curr_extra_vals, curr_causal_idxs, curr_bools)
                        if unmatched_idx < safe_index(val, indices, (indices == causal_idx) & ~curr_bool).shape[0]
                    ])
                    subgrouped_extra.append(np.nanmean(curr_extra_unmatched))

                if curr_extra_vals_2 is not None:
                    curr_extra_unmatched_2 = np.array([
                        safe_index(val, indices, (indices == causal_idx) & ~curr_bool)[unmatched_idx]
                        for val, indices, curr_bool in zip(curr_extra_vals_2, curr_causal_idxs, curr_bools)
                        if unmatched_idx < safe_index(val, indices, (indices == causal_idx) & ~curr_bool).shape[0]
                    ])
                    subgrouped_extra_2.append(np.nanmean(curr_extra_unmatched_2))

        grouped_vals.append(np.array(subgrouped_vals, dtype=object))
        grouped_sem.append(np.array(subgrouped_sem, dtype=object))
        grouped_bools.append(np.array(subgrouped_bools, dtype=object))
        if curr_metrics is not None:
            grouped_metrics.append(np.array(subgrouped_metrics, dtype=object))
            grouped_metrics_sem.append(np.array(subgrouped_metrics_sem, dtype=object))
        if curr_extra_vals is not None:
            grouped_extra.append(np.array(subgrouped_extra, dtype=object))
        if curr_extra_vals_2 is not None:
            grouped_extra_2.append(np.array(subgrouped_extra_2, dtype=object))

    if curr_metrics is not None:
        return grouped_params, grouped_vals, grouped_sem, grouped_bools, grouped_metrics, grouped_metrics_sem, grouped_extra, grouped_extra_2
    return grouped_params, grouped_vals, grouped_sem, grouped_bools, grouped_extra, grouped_extra_2


@click.command()
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("observed_hap", type=click.Path(path_type=Path))
@click.argument("causal_hap", type=click.Path(path_type=Path))
@click.option(
    "-f",
    "--files",
    type=click.Path(path_type=Path),
    default=None,
    show_default=True,
    help="A file listing a subset of observed hap files to which we will apply the match",
)
@click.option(
    "-m",
    "--metrics",
    type=Path,
    default=None,
    show_default=True,
    help="Plot fine-mapping metrics, as well",
)
@click.option(
    "--use-metric",
    type=str,
    default=None,
    show_default="observed LD",
    help="Use this fine-mapping metric to color the plot",
)
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all genotypes",
    help="""
    The region from which to extract genotypes; ex: 'chr1:1234-34566' or 'chr7'\n
    For this to work, the VCF must be indexed and the seqname must match!""",
)
@click.option(
    "--maf",
    type=float,
    default=None,
    show_default="no filtering",
    help="Ignore variants with a MAF below this threshold",
)
@click.option(
    "-i",
    "--observed-id",
    type=str,
    show_default="the best haplotype in the file",
    help="A haplotype ID from the .hap file output by happler",
)
@click.option(
    "-c",
    "--causal-id",
    type=str,
    show_default="the first haplotype in the file",
    help="A haplotype ID from the causal .hap file",
)
@click.option(
    "--hide-extras",
    is_flag=True,
    default=False,
    show_default=True,
    help="Hide the observed haps that have no causal hap match",
)
@click.option(
    "--pickle-out",
    is_flag=True,
    default=False,
    show_default=True,
    help="Save the output as a pickle file as well",
)
@click.option(
    "--order",
    type=str,
    default=None,
    show_default="The order of the wildcards in the observed_hap path",
    help="The order of the parameters in the plot, as a comma-separated list",
)
@click.option(
    "--phenos",
    type=Path,
    default=None,
    show_default=True,
    help="If passed, the optimal BF threshold is also plotted",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A PNG file containing the desired heatmap",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="ERROR",
    show_default=True,
    help="The level of verbosity desired",
)
def main(
    genotypes: Path,
    observed_hap: Path,
    causal_hap: Path,
    files: Path = None,
    metrics: Path = None,
    use_metric: str = None,
    region: str = None,
    maf: float = None,
    observed_id: str = None,
    causal_id: str = None,
    hide_extras: bool = False,
    pickle_out: bool = False,
    order: str = None,
    phenos: bool = False,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "ERROR",
):
    """
    Create a plot summarizing the LD between the observed and causal haplotypes across
    a range of parameters

    Parameters are inferred from brace expressions injected into the paths to the files
    For example, a path like "out/19_45401409-46401409/happler/hap/{beta}/happler.hap"
    will infer the beta values from the portion of the path labeled '{beta}'
    """
    log = logging.getLogger("parameter_plot", level=verbosity)

    # initialize the genotypes class
    gts = GenotypesVCF
    if genotypes.suffix == ".pgen":
        gts = GenotypesPLINK
    gts = gts(genotypes, log=log)

    # which files should we originally consider?
    # by default, we just grab as many as we can
    if files is not None:
        log.debug("Obtaining filtering list")
        with open(files, "r") as files_subset_file:
            files = files_subset_file.read().splitlines()

    # extract parameters and parameter values by globbing wildcards
    # params will be a dictionary mapping parameter names to lists of values
    params = dict(glob_wildcards(observed_hap, files=files)._asdict())
    # convert the dictionary to a numpy mixed dtype array
    dtypes = {k: DTYPES[k] for k in params.keys()}
    # if the user specified an order, then we will sort the parameters by that order
    if order is not None:
        new_params = {}
        for k in order.split(",")[::-1]:
            if k not in dtypes:
                log.warning(f"Parameter {k} not found in the observed hap path")
                continue
            dtypes = {k: dtypes.pop(k), **dtypes}
        for k in dtypes:
            new_params[k] = params.pop(k)
        params = new_params
    params = np.array(list(zip(*params.values())), dtype=list(dtypes.items()))
    params.sort()

    get_hap_fname = lambda hap_path, param_set: Path(str(hap_path).format(**dict(zip(dtypes.keys(), param_set))))

    # compute LD between the causal hap and the best observed hap across param vals
    all_ld_values = zip(*tuple(
        get_best_ld(
            gts,
            get_hap_fname(observed_hap, params[idx]),
            get_hap_fname(causal_hap, params[idx]),
            region=region,
            maf=maf,
            observed_id=observed_id,
            causal_id=causal_id,
            pheno=get_hap_fname(phenos, params[idx]),
            log=log
        )
        for idx in range(len(params))
    ))
    optimize_thresh = phenos is not None
    if optimize_thresh:
        ld_vals, ld_extras_idxs, ld_extras_bool, best_bf, causal_best_bf = all_ld_values
        best_bf = np.array(best_bf, dtype=object)
        causal_best_bf = np.array(causal_best_bf, dtype=object)
    else:
        ld_vals, ld_extras_idxs, ld_extras_bool = all_ld_values
        best_bf = None
    ld_vals = np.array(ld_vals, dtype=object)
    ld_extras_idxs = np.array(ld_extras_idxs, dtype=object)
    ld_extras_bool = np.array(ld_extras_bool, dtype=object)

    with open(output.with_suffix(".causal.pickle"), "wb") as f:
        pickle.dump([causal_best_bf], f)
    

    # extract fine-mapping metrics for the observed hap
    if metrics is not None:
        metrics = np.array([
            get_finemap_metrics(
                get_hap_fname(metrics, params[idx]),
                log=log
            )
            for idx in range(len(params))
        ], dtype=object)
        params, ld_vals, ld_sem, ld_extras_bool, metrics_vals, metrics_vals_sem, best_bf_grouped, causal_best_bf_grouped = group_by_rep(
            params, ld_vals, ld_extras_idxs, ld_extras_bool, metrics, extra_vals=best_bf, extra_vals_2=causal_best_bf,
        )
        extract_metrics_vals = lambda x, idx: np.array([i[:, idx] for i in x], dtype=object)
        metrics = {
            name: (
                extract_metrics_vals(metrics_vals, idx),
                extract_metrics_vals(metrics_vals_sem, idx)
            ) for idx, name in enumerate(metrics[0].dtype.names)
        }
    else:
        params, ld_vals, ld_sem, ld_extras_bool, best_bf_grouped, causal_best_bf_grouped = group_by_rep(
            params, ld_vals, ld_extras_idxs, ld_extras_bool, extra_vals=best_bf, extra_vals_2=causal_best_bf,
        )
    del dtypes["rep"]

    if pickle_out:
        with open(output.with_suffix(".pickle"), "wb") as f:
            if metrics is not None:
                pickle.dump([params, ld_vals, ld_sem, ld_extras_bool, metrics, best_bf_grouped, causal_best_bf_grouped], f)
            else:
                pickle.dump([params, ld_vals, ld_sem, ld_extras_bool, best_bf_grouped, causal_best_bf_grouped], f)

    diff_dtype = remove_same_valued_columns(params).dtype
    if use_metric is None:
        if len(diff_dtype) > 1:
            fig = plot_params(
                params,
                ld_vals,
                ld_sem,
                "Observed LD",
                ld_extras_bool,
                metrics=metrics,
                hide_extras=hide_extras,
                best_bf_vals=best_bf_grouped,
            )
        elif len(diff_dtype) == 1:
            fig = plot_params_simple(
                params,
                ld_vals,
                ld_sem,
                "Observed LD",
                ld_extras_bool,
                hide_extras=hide_extras,
            )
    else:
        if len(diff_dtype) > 2:
            fig = plot_params(
                params,
                metrics[use_metric][0],
                metrics[use_metric][1],
                use_metric,
                ld_extras_bool,
                hide_extras=hide_extras,
            )
        elif len(diff_dtype) == 1:
            fig = plot_params_simple(
                params,
                metrics[use_metric][0],
                metrics[use_metric][1],
                use_metric,
                ld_extras_bool,
                ld_extras_bool,
                hide_extras=hide_extras,
            )

    fig.tight_layout()
    fig.savefig(output,
    bbox_inches="tight", pad_inches=0.05)

if __name__ == "__main__":
    main()
