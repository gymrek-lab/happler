#!/usr/bin/env python
import pickle
from pathlib import Path
from logging import Logger

import click
import matplotlib
import numpy as np
matplotlib.use('Agg')
import numpy.typing as npt
import matplotlib.pyplot as plt

from haptools.logging import getLogger

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
LEGENDSIZE = 10
LABELSIZE = 12

# COMMAND-LINE:
# workflow/scripts/final-plots.py out/19_45401409-46401409/{sampsize}samples/plots/ld_range/{happler_params,finemapping_metrics}.pickle figures


def PlotBeta(datavals, colors, ldvals, labels, ax):
    for j in range(len(datavals)):
        data = datavals[j]
        idx = [i for i in range(len(data[0])) if (data[0][i]["ld"]==ldvals[j] and data[0][i]["num_haps"]==1)]
        betas = []
        means = []
        sds = []
        for i in idx:
            if len(data[1][i]) > 0:
                betas.append(data[0][i]["beta"])
                means.append(data[1][i][0])
                sds.append(data[2][i][0])
        ax.errorbar(betas, means, yerr=sds, color=colors[j], marker="o", ls='none', label=labels[j])
    ax.legend(loc="lower right", fontsize=LEGENDSIZE)
    ax.set_xlabel("Beta", size=LABELSIZE)
    ax.set_ylabel("Observed LD", size=LABELSIZE)

def PlotLD(datavals, colors, betavals, labels, ax):
    for j in range(len(datavals)):
        data = datavals[j]
        idx = [i for i in range(len(data[0])) if (data[0][i]["beta"]==betavals[j] and data[0][i]["num_haps"]==1)]
        ldvals = []
        means = []
        sds = []
        for i in idx:
            if len(data[1][i]) > 0:
                ldvals.append(data[0][i]["ld"])
                means.append(data[1][i][0])
                sds.append(data[2][i][0])
        ax.errorbar(ldvals, means, yerr=sds, color=colors[j], marker="o", ls='none', label=labels[j])
    ax.legend(loc="lower right", fontsize=LEGENDSIZE)
    ax.set_xlabel("Causal LD", size=LABELSIZE)
    ax.set_ylabel("Observed LD", size=LABELSIZE)

def PlotNumHaps(datavals, colors, ldvals, betavals, labels, ax):
    for j in range(len(datavals)):
        data = datavals[j]
        added_label = False
        for hapnum in range(1, 5):
            try:
                idx = np.where((data[0]['beta'] == betavals[j]) & (data[0]['num_haps'] == hapnum) & (data[0]['ld'] == ldvals[j]))[0][0]
            except IndexError as e:
                if ldvals[j] in np.unique(data[0]['ld']).tolist():
                    raise e
                else:
                    raise ValueError("Check that the desired LD vals exist in the data")
            for hn in range(3):
                if hn >= len(data[1][idx]): continue
                mean = data[1][idx][hn]
                sd = data[2][idx][hn]
                passed = data[3][idx][hn]
                if passed:
                    alpha = 1
                else: alpha = 0.25
                if not added_label:
                    ax.errorbar([hapnum+hn*0.1], [mean], yerr=[sd], color=colors[j], marker="o", ls='none', alpha=alpha, label=labels[j])
                    added_label = True
                else:
                    ax.errorbar([hapnum+hn*0.1], [mean], yerr=[sd], color=colors[j], marker="o", ls='none', alpha=alpha)
    ax.legend(loc="lower right", fontsize=LEGENDSIZE)
    ax.set_ylabel("Observed LD", size=LABELSIZE)
    ax.set_xlabel("# causal haplotypes", size=LABELSIZE)

def PlotAllPIPs(data, colors, ldvals, beta, labels, ax, include="all", alpha=1, legend=True):
    for j in range(len(colors)):
        added_label = False
        for hapnum in range(1, 4):
            idx = [i for i in range(len(data[0])) if (data[0][i]["beta"]==beta and \
                                                      data[0][i]["num_haps"]==hapnum and \
                                                      data[0][i]["ld"]==ldvals[j])][0]
            for hn in range(3):
                if hn >= len(data[1][idx]): continue
                mean = data[4]["pip"][0][idx][hn]
                sd = data[4]["pip"][1][idx][hn]
                passed = data[3][idx][hn]
                if passed:
                    alpha = 1
                else: alpha = 0.25
                if not added_label:
                    ax.errorbar([hapnum+hn*0.1], [mean], yerr=[sd], color=colors[j], marker="o", ls='none', alpha=alpha, label=labels[j])
                    added_label = True
                else:
                    ax.errorbar([hapnum+hn*0.1], [mean], yerr=[sd], color=colors[j], marker="o", ls='none', alpha=alpha)
    ax.set_xticks(range(1, 4))
    ax.set_xticklabels(range(1, 4))
    ax.set_ylabel("PIP", size=LABELSIZE)
    ax.set_xlabel("# causal haplotypes", size=LABELSIZE);

def PlotPIPSNoMatchVsMatch(data, colors, ldvals, beta, labels, ax):
    for j in range(len(colors)):
        pips_match = []
        pips_nomatch = []
        for hapnum in range(1, 4):
            idx = [i for i in range(len(data[0])) if (data[0][i]["beta"]==beta and \
                                                      data[0][i]["num_haps"]==hapnum and \
                                                      data[0][i]["ld"]==ldvals[j])][0]
            for hn in range(3):
                if hn >= len(data[1][idx]): continue
                mean = data[4]["pip"][0][idx][hn]
                sd = data[4]["pip"][1][idx][hn]
                passed = data[3][idx][hn]
                if passed:
                    pips_match.append(mean)
                else:
                    pips_nomatch.append(mean)
        ax.hist(pips_match, color=colors[j], alpha=1)
        ax.hist(pips_nomatch, color=colors[j], alpha=0.05)
    ax.set_xlabel("PIP", size=LABELSIZE)
    ax.set_ylabel("Frequency", size=LABELSIZE);

def PlotBestPipByBetaLDVals(data, colors, ldvals, labels, ax, include="all", alpha=1, legend=True):
    for j in range(len(colors)):
        idx = [i for i in range(len(data[0])) if (data[0][i]["ld"]==ldvals[j] and data[0][i]["num_haps"]==1)]
        betas = []
        means = []
        sds = []
        for i in idx:
            if len(data[1][i]) > 0:
                pips = list(data[4]["pip"][0][i])
                best_pip_idx = pips.index(max(pips))
                passed = data[3][i][best_pip_idx]
                if include == "nomatch" and passed:
                    continue
                if include == "match" and not passed:
                    continue
                betas.append(data[0][i]["beta"])
                means.append(data[4]["pip"][0][i][best_pip_idx])
                sds.append(data[4]["pip"][1][i][best_pip_idx])
        ax.errorbar(betas, means, yerr=sds, color=colors[j], marker="o", label=labels[j], alpha=alpha)
    if legend:
        ax.legend(loc="upper left", fontsize=LEGENDSIZE)
        ax.set_xlabel("Beta", size=LABELSIZE)
        ax.set_ylabel("Best haplotype PIP", size=LABELSIZE);
        #ax.set_title("1 causal haplotype")

def Fig2(data_3200, data_800, output_dir):
    LOWN_COLOR = "darkblue"
    HIGHN_COLOR = "darkred"
    LOWN_LDVAL = 0.73
    HIGHN_LDVAL = 0.71
    LOWN_BETA = 0.15
    HIGHN_BETA = 0.15

    fig, axs = plt.subplot_mosaic("ABC", figsize=(12,4))
    for n, (key, ax) in enumerate(axs.items()):
        if n == 0:
            PlotBeta([data_800, data_3200], [LOWN_COLOR, HIGHN_COLOR], [LOWN_LDVAL, HIGHN_LDVAL],
                    ["n=800; LD=%.2f"%LOWN_LDVAL, "n=3200; LD=%.2f"%HIGHN_LDVAL], ax)
        if n == 1:
            PlotLD([data_800, data_3200], [LOWN_COLOR, HIGHN_COLOR], [LOWN_BETA, HIGHN_BETA],
                ["n=800; B=%.2f"%LOWN_BETA, "n=3200; B=%.2f"%HIGHN_BETA], ax)
        if n == 2:
            PlotNumHaps([data_800, data_3200], [LOWN_COLOR, HIGHN_COLOR], [LOWN_LDVAL, HIGHN_LDVAL], [LOWN_BETA, HIGHN_BETA],
                        ["n=800; LD=%.2f; B=%.2f"%(LOWN_LDVAL, LOWN_BETA), "n=3200; LD=%.2f; B=%.2f"%(HIGHN_LDVAL, HIGHN_BETA)], ax)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        # Add figure subplot letters
        ax.text(-0.1, 1.1, key, transform=ax.transAxes,
                size=20, weight='bold')
    fig.tight_layout()
    fig.savefig(output_dir / "Happler_Fig2A_Draft.pdf")

def Fig3(data_3200, output_dir):
    LDVALS = [0.02, 0.26, 0.71, 0.94]
    colors = ["red", "orange", "blue", "purple"]
    BETA = 0.40
    fig, axs = plt.subplot_mosaic("ABC", figsize=(12,4))
    assert set(LDVALS) < set(np.unique(data_3200[0]["ld"]).tolist()), "Check that the desired LD vals exist in the data"
    for n, (key, ax) in enumerate(axs.items()):
        if n == 0:
            PlotBestPipByBetaLDVals(data_3200, colors, LDVALS, ["LD=%.2f"%item for item in LDVALS], ax=ax, include="all")
        #if n == 1:
        #    PlotAllPIPs(data_3200, colors, LDVALS, BETA, ["LD=%.2f"%item for item in LDVALS], ax=ax)
        #if n == 1:
        #    PlotPIPSNoMatchVsMatch(data_3200, colors, LDVALS, BETA, ["LD=%.2f"%item for item in LDVALS], ax=ax)
        if n == 1 or n == 2:
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        # Add figure subplot letters
        ax.text(-0.1, 1.1, key, transform=ax.transAxes,
                size=20, weight='bold')
    fig.tight_layout()
    fig.savefig(output_dir / "Happler_Fig3_Draft.pdf")


@click.command()
@click.argument("sim_params", type=click.Path(path_type=Path))
@click.argument("sim_params_finemap", type=click.Path(path_type=Path))
@click.argument("output_dir", type=click.Path(path_type=Path))
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="DEBUG",
    show_default=True,
    help="The level of verbosity desired",
)
def main(
    sim_params: Path,
    sim_params_finemap: Path,
    output_dir: Path,
    verbosity: str = "DEBUG",
):
    """
    Make plots for the happler paper
    """
    log = getLogger("final-plots", verbosity)

    output_dir.mkdir(parents=True, exist_ok=True)

    with open(str(sim_params).format(sampsize=800), "rb") as f:
        data_800 = pickle.load(f)
        ldvals = np.unique(data_800[0]["ld"]).tolist()
        log.info(f"LD for 800: {ldvals}")
    with open(str(sim_params).format(sampsize=3200), "rb") as f:
        data_3200 = pickle.load(f)
        ldvals = np.unique(data_3200[0]["ld"]).tolist()
        log.info(f"LD for 3200: {ldvals}")
    Fig2(data_3200, data_800, output_dir)
    with open(str(sim_params_finemap).format(sampsize=3200), "rb") as f:
        data_3200 = pickle.load(f)
    Fig3(data_3200, output_dir)


if __name__ == "__main__":
    main()
