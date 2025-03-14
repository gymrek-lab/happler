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
# this gets imported later down only if --pos-type is specified
# from sklearn.metrics import roc_curve, auc

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
        converters={"P": lambda val: Decimal(val if val != "NA" else 1)},
    ).rename(columns=PLINK_COLS)
    df = df.sort_values("pos")
    df["pval"] = df["pval"].fillna(np.inf)
    return df


def get_pval(
    linear: Path,
    snp_id: str,
    log: logging.Logger = None
) -> float:
    """
    Extract the p-value from the linear file which corresponds with the last SNP ID in
    the snplist file

    Parameters
    ----------
    linear: Path
        The path to a PGEN file containing genotypes for all haplotypes and their SNPs
    snp_id: str
        The ID of the target SNP
    log: logging.Logger, optional
        A logging object to write any debugging and error messages

    Returns
    -------
    float
        The -log10 p-value (from the linear file) of the last SNP in the snplist file
    """
    # load the linear file
    df = load_linear_file(linear)
    pval = df[df.id == snp_id].iloc[0]["pval"]
    pval = np.float64(-np.log10(pval))
    return pval


def get_snp_id(
    snplist: Path,
    log: logging.Logger = None
) -> float:
    """
    Extract the last SNP ID from the snplist file

    Parameters
    ----------
    snplist: Path
        The path to a .hap file containing a set of haplotypes
    log: logging.Logger, optional
        A logging object to write any debugging and error messages

    Returns
    -------
    str
        The SNP ID from the snplist file
    """
    # load the snp ID from the snplist file
    if snplist.suffix == ".snplist":
        with open(snplist) as snplist_file:
            snp_id = snplist_file.readlines()[-1].split("\t")[0]
    else:
        # load from a hap file
        with open(snplist) as hap_file:
            snp_id = [
                line.split("\t")[4]
                for line in hap_file.readlines()
                if line.startswith("V")
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
    xmax = np.max(np.abs(x))
    ymax = np.max(np.abs(y))
    xbinwidth = (xmax/25)/4
    ybinwidth = (ymax/25)/4
    padding = 50
    if zoom:
        ax.set_xlim(0-(xmax/padding), xmax+(xmax/padding))
        ax.set_ylim(0-(ymax/padding), ymax+(ymax/padding))
    xlim = (int(xmax/xbinwidth) + 1) * xbinwidth
    ylim = (int(ymax/ybinwidth) + 1) * ybinwidth
    xbins = np.arange(0, xlim + xbinwidth, xbinwidth)
    ybins = np.arange(0, ylim + ybinwidth, ybinwidth)
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
    help="Values originating from files with this wildcard will be colored the same",
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
    help="The significance threshold"
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
    {region}.{type}/1/out.{rep}.{name}.glm.linear
    {region}.{type}/2/out.{rep}.{name}.glm.linear

    The snplists files can either be formatted as snplist files or hap files
    """
    log = getLogger("midway_manhattan_summary", level=verbosity)

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
    if color is not None:
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

    # color reps the same
    if color is not None:
        log.debug("Setting up colors")
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
        if color is not None:
            colors = np.delete(colors, exclude_idxs, 0)

    # check that everything is still kosher wrt to the rep field
    if 'rep' in dtypes:
        assert (params[axes_idxs[ax_labs[0]]]['rep'] == params[axes_idxs[ax_labs[1]]]['rep']).all()

    # compute -log10 pval for the target SNP in each linear file
    # Note: this will return a 2D array where each column corresponds to a plot axis
    log.debug(f"Extracting pvals from linear files")
    vals = np.array(
        [
            [
                get_pval(
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

    # remove any vals that were NA (but got converted to 0)
    # and also any vals that were greater than max-val
    na_rows = (vals == 0).any(axis=1) | (vals > max_val).any(axis=1)
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
        log.debug("Plotting AUROC")
        from sklearn.metrics import roc_curve, auc
        # first, get the roc curve vals
        pos_labs = np.zeros(vals.shape, dtype=np.bool_)
        pos_labs[:,pos_type] = True
        fpr, tpr, threshold = roc_curve(pos_labs.flatten("F"), vals.flatten("F"), drop_intermediate=True)
        if thresh is None:
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
            thresh = 10**(-threshold[thresh_idx])
            print(f"Threshold: {thresh}")
        else:
            thresh_idx = np.argmax(threshold < -np.log10(thresh))
        roc_auc = auc(fpr, tpr)
        print(f"AUROC: {roc_auc}")
        # now, make the fig
        fig = plt.figure(figsize=(11, 6), layout='constrained')
        subfigs = fig.subfigures(1, 2, wspace=0, width_ratios=(6, 5))
        subfig = subfigs[0]
        # make the roc plot
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
    # Draw the scatter plot and marginals.
    log.debug("Creating scatter plot")
    if color is None:
        scatter_hist(vals[:,1], vals[:,0], ax, ax_histx, ax_histy)
    else:
        scatter_hist(vals[:,1], vals[:,0], ax, ax_histx, ax_histy, colors=colors)
    max_val = vals.max()
    fig.text(0.98, 0.98, f'Threshold: {thresh:.2f}', ha='right', va='top', fontsize=15)
    thresh = -np.log10(thresh)
    ax.set_xlabel(case_type + ": " + ax_labs[1])
    ax.set_ylabel(case_type + ": " + ax_labs[0])
    ax.axline((0,0), (max_val, max_val), linestyle="--", color="orange")
    ax.axline((0,thresh), (thresh, thresh), color="red")
    ax_histx.axline((thresh,0), (thresh, thresh), color="red")
    ax.axline((thresh,0), (thresh, thresh), color="red")
    ax_histy.axline((0,thresh), (thresh, thresh), color="red")
    ax_histy.spines['top'].set_visible(False)
    ax_histx.spines['top'].set_visible(False)
    ax_histy.spines['right'].set_visible(False)
    ax_histx.spines['right'].set_visible(False)
    plt.savefig(output)


if __name__ == "__main__":
    main()
