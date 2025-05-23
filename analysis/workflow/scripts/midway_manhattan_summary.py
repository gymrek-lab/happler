#!/usr/bin/env python
import math
import logging
import warnings
from pathlib import Path
from decimal import Decimal
from functools import partial

import click
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# this gets imported later down only if --pos-type is specified
# from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score, precision_recall_fscore_support

from snakemake_io import glob_wildcards, remove_regexes


def getLogger(name: str = None, level: str = "ERROR", exact_time=False):
    """
    Retrieve a Logger object

    Parameters
    ----------
    name : str, optional
        The name of the logging object
    level : str, optional
        The level of verbosity for the logger
    """
    if name is None:
        name = ""
    else:
        name = "." + name

    # create logger
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(level)

    # create formatter
    db_time = (
        ("|%(asctime)s" + (".%(msecs)03d" if exact_time else ""))
        if level == "DEBUG"
        else ""
    )
    formatter = logging.Formatter(
        fmt="[%(levelname)8s" + db_time + "] %(message)s (%(filename)s:%(lineno)s)",
        datefmt="%H:%M:%S",
    )

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)

    return logger


PLINK_COLS = {
    "#CHROM": "chromosome",
    "POS": "pos",
    "ID": "id",
    "P": "pval",
}


def load_linear_file(linear_fname: Path):
    keep_cols = list(PLINK_COLS.keys())
    df = pd.read_csv(
        linear_fname,
        sep="\t",
        header=0,
        usecols=keep_cols,
        converters={"P": lambda val: Decimal(val) if val != "NA" else np.nan},
    ).rename(columns=PLINK_COLS)
    # there should only be NA values if VIF is infinite or something like that
    return df


def get_finemap_metrics(
    metrics_path: Path,
    keep_hap_ids: bool = False,
    log: logging.Logger = None
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
    # copy of this function in parameter_plot.py 
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
        metrics = np.delete(metrics, 0, axis=1)
    return metrics


def get_pval(
    linear: Path,
    snp_id: str,
    bic: bool = False,
    log: logging.Logger = None
) -> float:
    """
    Extract the p-value from the linear file which corresponds with the last SNP ID in
    the snplist file

    Parameters
    ----------
    linear: Path
        The path to a .linear file containing p-values for a bunch of SNPs
    snp_id: str
        The ID of the target SNP
    bic: bool, optional
        Whether to treat the pval as a delta BIC value, instead
    log: logging.Logger, optional
        A logging object to write any debugging and error messages

    Returns
    -------
    float
        The -log10 p-value (from the linear file) of the target SNP
    """
    # load the linear file
    df = load_linear_file(linear)
    pval = df[df.id == snp_id].iloc[0]["pval"]
    if bic:
        pval = math.log(pval)
    else:
        pval = -np.log10(pval)
    return np.float64(pval)


def get_pip(
    metrics_file: Path,
    hap_id: str,
    target_metric: str = "pip",
    log: logging.Logger = None,
) -> float:
    """
    Extract the PIP from the finemapping metrics file which corresponds with the hap ID

    Parameters
    ----------
    metrics_file: Path
        The path to a susie_metrics.tsv file containing fine-mapping metrics
    hap_id: str
        The ID of the target hap
    target_metric: str, optional
        The name of the metric that we want to load
    log: logging.Logger, optional
        A logging object to write any debugging and error messages

    Returns
    -------
    float
        The PIP of the target hap
    """
    # load the metrics file
    metrics = get_finemap_metrics(metrics_file, keep_hap_ids=True, log=log)
    # check: were there any credible sets?
    if len(metrics) == 1 and not metrics["num_credible_sets"][0]:
        # if there weren't any, just output the default value
        return metrics[target_metric][0]
    else:
        pip = metrics[["hap_id", target_metric]]
        hap_idxs = pip["hap_id"] == hap_id
        # could we find the hap in the metrics file?
        # if we can't find the hap, it's because it wasn't included in the input dataset
        if (sum(hap_idxs) != 0):
            pip = pip[hap_idxs]
        # there should only be one value but we use np.mean in case it couldn't find the hap
        return np.mean(pip[target_metric])


def get_snp_id(
    snplist: Path,
    get_causal: bool = False,
    hap_id: bool = False,
    log: logging.Logger = None
) -> float:
    """
    Extract the last SNP ID from the snplist file

    Parameters
    ----------
    snplist: Path
        The path to a .hap file containing a set of haplotypes
    get_causal: bool, optional
        If True, also check whether the file has "# not.causal" in it
    hap_id: bool, optional
        If True and the file is a hap file, return the hap ID instead of the variant ID
    log: logging.Logger, optional
        A logging object to write any debugging and error messages

    Returns
    -------
    str|bool
        The SNP ID from the snplist file if get_causal is False. Otherwise, a bool
        indicating whether the file has "# not.causal" in it
    """
    # load the snp ID from the snplist file
    if snplist.suffix == ".snplist":
        with open(snplist) as snplist_file:
            if get_causal:
                snp_id = ("# not.causal" in snplist_file.read().splitlines())
            else:
                snp_id = snplist_file.readlines()[-1].split("\t")[0]
    else:
        # load from a hap file
        with open(snplist) as hap_file:
            if get_causal:
                snp_id = ("# not.causal" in hap_file.read().splitlines())
            else:
                snp_id = [
                    line.split("\t")[4]
                    for line in hap_file.readlines()
                    if line.startswith("H" if hap_id else "V")
                ][-1]
    return snp_id


def scatter_hist(x, y, ax, ax_histx, ax_histy, colors=None, zoom=True):
    """
    Adapted from https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html#sphx-glr-gallery-lines-bars-and-markers-scatter-hist-py
    """
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    if colors is None:
        ax.scatter(x, y)
    else:
        ax.scatter(x, y, color=colors)

    # now determine nice limits by hand:
    xmin = np.min(x)
    ymin = np.min(y)
    xmax = np.max(x)
    ymax = np.max(y)
    xbinwidth = ((xmax-xmin)/25)/2
    ybinwidth = ((ymax-ymin)/25)/2
    if not xbinwidth:
        xbinwidth = 1
    if not ybinwidth:
        ybinwidth = 1
    padding = 50
    if zoom:
        ax.set_xlim(xmin-(xmax/padding), xmax+(xmax/padding))
        ax.set_ylim(ymin-(ymax/padding), ymax+(ymax/padding))
    xlim = (int(xmax/xbinwidth) + 1) * xbinwidth
    ylim = (int(ymax/ybinwidth) + 1) * ybinwidth
    xbins = np.arange(xmin, xlim + xbinwidth, xbinwidth)
    ybins = np.arange(ymin, ylim + ybinwidth, ybinwidth)
    ax_histx.hist(x, bins=xbins)
    ax_histy.hist(y, bins=ybins, orientation='horizontal')


@click.command()
@click.argument("linears", type=click.Path(path_type=Path))
@click.argument("snplists", type=click.Path(path_type=Path))
@click.argument("case_type", type=str)
@click.option(
    "-f",
    "--files",
    type=click.Path(path_type=Path),
    default=None,
    show_default=True,
    help="A file containing a subset of linear files to which we will apply the linears match",
)
@click.option(
    "-e",
    "--exclude",
    type=click.Path(path_type=Path),
    default=None,
    show_default=True,
    help="A file containing a list of linear files to exclude from the analysis",
)
@click.option(
    "--max-val",
    type=float,
    default=np.inf,
    show_default=True,
    help="Don't show p-values larger than this -log10 val"
)
@click.option(
    "-c",
    "--color",
    type=str,
    default=None,
    show_default=True,
    help=(
        "Values originating from files with this wildcard will be colored the same. Or"
        " if 'not.causal', then snplist/hap files with '# not.causal' in them will be colored."
    ),
)
@click.option(
    "-p",
    "--pos-type",
    type=str,
    default=None,
    show_default="",
    help="If you want an ROC plot, provide the case-type that represents positives",
)
@click.option(
    "-t",
    "--thresh",
    type=float,
    default=None,
    show_default="FDR 0.05",
    help="The significance threshold. Inferred at an FDR of 0.05 if not specified"
)
@click.option(
    "--kind",
    type=click.Choice(["pval", "bic", "pip", "cs_length", "in_credible_set", "num_credible_sets", "cs_purity", "has_highest_pip", "best_variant_pip"]),
    default="pval",
    show_default=True,
    help="What kind of value are we plotting?",
)
@click.option(
    "--no-log10",
    is_flag=True,
    show_default=True,
    default=False,
    help="Do not show points on the -log10 scale but on the normal scale instead",
)
@click.option(
    "--reverse",
    is_flag=True,
    show_default=True,
    default=False,
    help="If True, put a negative sign in front of each (p)val. This flips the ROC and PRC upside down",
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
    linears: Path,
    snplists: Path,
    case_type: str,
    files: Path = None,
    exclude: Path = None,
    max_val: float = np.inf,
    color: str = None,
    pos_type: str = None,
    thresh: float = None,
    kind: str = "pval",
    no_log10: bool = False,
    reverse: bool = False,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Plot p-values from midway-manhattan simulations in two different situations (case
    - types) against each other

    Each file is inferred from brace expressions injected into the paths to the files
    For example, a path like "{region}.{type}/{tswitch}/out.{rep}.{name}.glm.linear"
    will infer the wildcards and match them between cases. Wildcards will be expanded
    across all inputs. So, for example, if the case_type is "tswitch", then these two
    files will be matched together:
    {region}.{type}/1/out.{rep}.{name}.linear
    {region}.{type}/2/out.{rep}.{name}.linear

    The snplists files can either be formatted as snplist files or hap files
    """
    log = getLogger("midway_manhattan_summary", level=verbosity)
    final_metrics={}

    # set the kind of value
    bic = False
    is_finemap_metric = False
    if kind == "bic":
        bic = True
    elif kind in ("pip", "cs_length", "in_credible_set", "num_credible_sets", "cs_purity", "has_highest_pip", "best_variant_pip"):
        is_finemap_metric = True

    # if this is a pval, take the -log10 of it, otherwise take the ln of it
    if kind == "pval":
        tsfm_pval = lambda pval: -np.log10(pval) if pval != 0 else np.inf
        rvrs_tsfm = lambda pval: np.power(10, -pval)
    elif bic:
        tsfm_pval = lambda val: np.log(val) if val != 0 else -np.inf
        rvrs_tsfm = lambda pval: np.power(10, pval)
    elif reverse:
        tsfm_pval = lambda val: -val
        rvrs_tsfm = lambda pval: -pval        
    else:
        tsfm_pval = lambda val: val
        rvrs_tsfm = lambda pval: pval

    # which files should we originally consider?
    # by default, we just grab as many as we can
    if files is not None:
        log.debug("Obtaining filtering list")
        with open(files, "r") as files_subset_file:
            files = files_subset_file.read().splitlines()

    # extract parameters and parameter values by globbing wildcards
    # params will be a dictionary mapping parameter names to lists of values
    params = dict(glob_wildcards(linears, files=files)._asdict())
    # check that case_type and color are wildcards
    assert case_type in params.keys()
    if color is not None and color != "not.causal":
        assert color in params.keys()
    dtypes = {k: "U30" for k in params.keys()}
    log.debug(f"Extracted paramter values {tuple(dtypes.keys())}")
    # convert the dictionary to a numpy mixed dtype array
    params = np.array(list(zip(*params.values())), dtype=list(dtypes.items()))
    params.sort()

    linears_wo_regexes = remove_regexes(str(linears))
    get_fname = lambda path, param_set: Path(str(path).format(**dict(zip(dtypes.keys(), param_set))))

    log.debug(f"Extracting SNP IDs from {len(params)} .snplist files")
    snp_IDs = {
        idx: get_snp_id(
            get_fname(str(snplists), params[idx]),
            hap_id=is_finemap_metric,
            log=log
        )
        for idx in range(len(params))
    }

    log.debug("Setting up axes")
    axes_idxs = {
        case_type_val: np.where(params[case_type] == case_type_val)[0]
        for case_type_val in tuple(np.unique(params[case_type]))
    }
    assert len(axes_idxs) == 2
    ax_labs = tuple(axes_idxs.keys())
    if pos_type is not None:
        assert pos_type in ax_labs
        pos_type = ax_labs.index(pos_type)

    # color by specific wildcard the same
    log.debug(f"Setting up colors variable for color: {color}")
    if color == "not.causal":
        causal_haps = {
            idx: get_snp_id(
                get_fname(str(snplists), params[idx]),
                get_causal=(color == "not.causal"),
                log=log
            )
            for idx in range(len(params))
        }
        color_vals = axes_idxs[ax_labs[0]]
        # orange = causal and blue = non-causal
        cmap = ("blue", "orange")
        colors = np.array([cmap[causal_haps[i]] for i in color_vals])
    elif color is not None:
        color_vals = params[axes_idxs[ax_labs[0]]][color]
        # check that the regions are the same in each
        assert (color_vals == params[axes_idxs[ax_labs[1]]][color]).all()
        color_vals_idxs = np.unique(color_vals)
        color_vals_idxs = dict(
            zip(color_vals_idxs, np.linspace(0, 1, len(color_vals_idxs)))
        )
        cmap = plt.get_cmap('turbo')
        colors = np.array(
            cmap([color_vals_idxs[reg] for reg in color_vals])
        )

    # remove any vals we were supposed to exclude
    if exclude is not None:
        with open(exclude, "r") as excludes_file:
            files_to_exclude = excludes_file.read().splitlines()
        exclude_params = dict(glob_wildcards(linears, files=files_to_exclude)._asdict())
        exclude_params = np.array(list(zip(*exclude_params.values())), dtype=list(dtypes.items()))
        assert len(exclude_params) == len(files_to_exclude)
        exclude_params.sort()
        exclude_idxs = np.searchsorted(params, exclude_params)
        assert (params[exclude_idxs] == exclude_params).all()
        exclude_idxs = np.concatenate(tuple(
            np.where(np.isin(axes_idxs[case_type_val], exclude_idxs))[0]
            for case_type_val in ax_labs
        ))
        log.debug(f"Excluding {len(exclude_idxs)} points")
        for ax_lab in ax_labs:
            axes_idxs[ax_lab] = np.delete(axes_idxs[ax_lab], exclude_idxs, 0)
        if color is not None and color != "not.causal":
            colors = np.delete(colors, exclude_idxs, 0)

    # check that everything is still kosher wrt to the rep field
    if 'rep' in dtypes:
        assert (params[axes_idxs[ax_labs[0]]]['rep'] == params[axes_idxs[ax_labs[1]]]['rep']).all()

    # compute -log10 pval for the target SNP in each linear file
    # Note: this will return a 2D array where each column corresponds to a plot axis
    log.debug(f"Extracting pvals from linear files")
    if is_finemap_metric:
        get_val_method = partial(get_pip, target_metric=kind)
    else:
        get_val_method = partial(get_pval, bic=bic)
    vals = np.array(
        [
            [
                get_val_method(
                    get_fname(linears_wo_regexes, params[idx]),
                    snp_IDs[idx],
                    log=log
                )
                for idx in axes_idxs[case_type_val]
            ]
            for case_type_val in axes_idxs
        ]
    ).T
    log.debug(f"Found {vals.shape[0]} points")

    # remove any rows that were nan or greater than max-val (which defaults to inf)
    if no_log10:
        if max_val == float("inf"):
            if not bic:
                max_val = 1
        na_rows = np.isnan(vals).any(axis=1) | (rvrs_tsfm(vals) >= max_val).any(axis=1)
    else:
        na_rows = np.isnan(vals).any(axis=1) | (vals >= max_val).any(axis=1)
    not_na_rows_idxs = {k: v[~na_rows] for k,v in axes_idxs.items()}
    if color is not None:
        colors = colors[~na_rows]
    not_na_rows_idxs = np.sort(np.concatenate(list(not_na_rows_idxs.values())))
    params = params[not_na_rows_idxs]
    vals = vals[~na_rows]
    log.debug(f"Removed {na_rows.sum()} NA points, resulting in {vals.shape[0]} remaining")

    # finally, make the plot
    # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
    # the size of the marginal Axes and the main Axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    # And if requested, compute an ROC curve, as well
    if pos_type is not None:
        log.debug("Plotting ROC and PRC")
        from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score, precision_recall_fscore_support
        # first, get the roc curve vals
        pos_labs = np.zeros(vals.shape, dtype=np.bool_)
        pos_labs[:,pos_type] = True
        y_true = pos_labs.flatten("F")
        y_score = vals.flatten("F")
        if reverse:
            y_score = -y_score
        fpr, tpr, roc_threshold = roc_curve(y_true, y_score, drop_intermediate=True)
        precision, recall, prc_threshold = precision_recall_curve(y_true, y_score, drop_intermediate=True)

        # now, try to get the optimal threshold
        with warnings.catch_warnings():
            # safe to ignore runtime warning caused by division of 0 by 0
            warnings.simplefilter("ignore")
            # compute FDR for each threshold
            # note that this computation only works because number of TN == number of TP
            fdr = fpr / (fpr + tpr)
            # set nan to 0
            fdr[np.isnan(fdr)] = 0
        # find the threshold (last index) where FDR <= 0.05
        thresh_idx = np.where(fdr <= 0.05)[0][-1]
        if bic:
            optimal_thresh = roc_threshold[thresh_idx]
        else:
            optimal_thresh = rvrs_tsfm(roc_threshold[thresh_idx])
        final_metrics["Significance Threshold"] = optimal_thresh
        if thresh is not None:
            thresh_idx = np.argmax(roc_threshold < tsfm_pval(thresh))
        else:
            thresh = optimal_thresh
        roc_auc = auc(fpr, tpr)
        prc_ap = average_precision_score(y_true, y_score)
        # Find the index where thresholds > log_thresh b/c prc_threshold increases from 0 to inf
        prc_thresh_idx = np.argmax(prc_threshold > tsfm_pval(thresh))

        # Ensure index is within bounds (prc_threshold is shorter with precision/recall than roc)
        if prc_thresh_idx >= len(precision):
            prc_thresh_idx = len(precision) - 1
        final_metrics["AUROC"] = roc_auc
        final_metrics["Average Precision"] = roc_auc

        # now, make the fig
        fig = plt.figure(figsize=(16, 6), layout='constrained')
        subfigs = fig.subfigures(1, 3, wspace=0, width_ratios=(6, 5, 5))
        # make the roc plot
        subfig = subfigs[0]
        ax_roc_gs = subfigs[1].add_gridspec(
            2, 1, height_ratios=(1, 5), wspace=0.01, hspace=0.01,
        )
        ax_roc = subfigs[1].add_subplot(ax_roc_gs[1, 0])
        ax_roc.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
        ax_roc.legend(loc = 'lower right')
        ax_roc.plot([0, 1], [0, 1], "--", color="orange")
        tpr_thresh = tpr[thresh_idx]
        if tpr_thresh != 0:
            ax_roc.axline((0, tpr_thresh), (tpr_thresh, tpr_thresh), color="red", lw=0.9)
        fpr_thresh = fpr[thresh_idx]
        if fpr_thresh != 0:
            ax_roc.axline((fpr_thresh, 0), (fpr_thresh, fpr_thresh), color="red", lw=0.9)
        ax_roc.set_xlim([-0.005, 1.005])
        ax_roc.set_ylim([-0.005, 1.005])
        ax_roc.set_ylabel('True Positive Rate')
        ax_roc.set_xlabel('False Positive Rate')
        final_metrics["FPR"] = fpr_thresh
        final_metrics["TPR"] = tpr_thresh
        # now, make the prc plot
        ax_prc_gs = subfigs[2].add_gridspec(
            2, 1, height_ratios=(1, 5), wspace=0.01, hspace=0.01,
        )
        ax_prc = subfigs[2].add_subplot(ax_prc_gs[1, 0])
        ax_prc.plot(recall, precision, 'b', label = 'AP = %0.2f' % prc_ap)
        ax_prc.legend(loc = 'lower right')
        ax_prc.plot([0, 1], [0, 1], "--", color="orange")
        precision_thresh = precision[prc_thresh_idx]
        if precision_thresh != 0:
            ax_prc.axline((0, precision_thresh), (precision_thresh, precision_thresh), color="red", lw=0.9)
        recall_thresh = recall[prc_thresh_idx]
        if recall_thresh != 0:
            ax_prc.axline((recall_thresh, 0), (recall_thresh, recall_thresh), color="red", lw=0.9)
        ax_prc.set_xlim([-0.005, 1.005])
        ax_prc.set_ylim([-0.005, 1.005])
        ax_prc.set_ylabel('Precision')
        ax_prc.set_xlabel('Recall')
        # compute final metrics
        precision, recall, fscore, support = precision_recall_fscore_support(y_true, y_score > tsfm_pval(thresh), pos_label=1)
        final_metrics["Precision"] = precision[1]
        final_metrics["Recall"] = recall[1]
        final_metrics["Fscore"] = fscore[1]
        final_metrics["Support"] = support[1]
    else:
        # just make the fig
        fig = plt.figure(figsize=(6, 6), layout='constrained')
        subfig = fig
    gs = subfig.add_gridspec(
        2, 2, width_ratios=(5, 1), height_ratios=(1, 5), wspace=0.01, hspace=0.01,
    )
    # Create the Axes.
    ax = subfig.add_subplot(gs[1, 0])
    ax_histx = subfig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = subfig.add_subplot(gs[1, 1], sharey=ax)
    # if --no-log10, then we should reverse the values back from the log scale
    if no_log10:
        vals = rvrs_tsfm(vals)
    # Draw the scatter plot and marginals.
    log.debug("Creating scatter plot")
    if color is None:
        scatter_hist(vals[:,1], vals[:,0], ax, ax_histx, ax_histy)
    else:
        scatter_hist(vals[:,1], vals[:,0], ax, ax_histx, ax_histy, colors=colors)
    if thresh is not None:
        threshold_type = "P-value"
        if bic:
            threshold_type = "Bayes factor"
        elif is_finemap_metric:
            if kind == "pip":
                threshold_type = "PIP"
            else:
                threshold_type = ""
        fig.text(0.98, 0.98, f'{threshold_type} threshold: {thresh:.2f}', ha='right', va='top', fontsize=15)
        if not no_log10 and not is_finemap_metric:
            thresh = tsfm_pval(thresh)
    ax.set_xlabel(case_type + ": " + ax_labs[1])
    ax.set_ylabel(case_type + ": " + ax_labs[0])
    ax.axline((0,0), (vals.max(), vals.max()), linestyle="--", color="orange")
    if thresh is not None and thresh != 0:
        ax.axline((0,thresh), (thresh, thresh), color="red")
        ax_histx.axline((thresh,0), (thresh, thresh), color="red")
        ax.axline((thresh,0), (thresh, thresh), color="red")
        ax_histy.axline((0,thresh), (thresh, thresh), color="red")
    ax_histy.spines['top'].set_visible(False)
    ax_histx.spines['top'].set_visible(False)
    ax_histy.spines['right'].set_visible(False)
    ax_histx.spines['right'].set_visible(False)
    plt.savefig(output)
    # print the final metrics to stdout
    print("\t".join(final_metrics.keys()))
    print("\t".join(map(str, final_metrics.values())))


if __name__ == "__main__":
    main()
