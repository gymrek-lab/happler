#!/usr/bin/env python
import re
import copy
from pathlib import Path
from logging import Logger

import click
import matplotlib
import numpy as np
matplotlib.use('Agg')
import numpy.typing as npt
import statsmodels.api as sm
import matplotlib.pyplot as plt

from haptools.logging import getLogger
from happler.tree import Haplotype as HapplerHaplotype
from happler.tree.assoc_test import AssocTestSimpleSM, AssocTestSimpleCovariates
from haptools.data import (
    Genotypes,
    Phenotypes,
    Haplotypes,
    GenotypesVCF,
    GenotypesPLINK,
)

FIGSIZE = 6

# TODO one day:
# Instead of using an orange dot to symbolize the haplotype, plot it as a black bar
# extending from its start to its end coordinates

def make_manhattan(
    ax: plt.Axes,
    positions: npt.NDArray,
    pvals: npt.NDArray,
    red_mask: npt.NDArray = None,
    exclude_mask: npt.NDArray = None,
    orange_mask: npt.NDArray = None,
):
    """
    Add a Manhattan plot to the desired Axes object

    Parameters
    ----------
    ax: plt.Axes
        The matplotlib axes to plot onto
    positions: npt.NDArray
        A 1D array of base pair positions for each SNP
    pvals: npt.NDArray
        A 1D array of p-values of the same length as positions
    red_mask: npt.NDArray, optional
        If provided, a bool array of pvals to highlight in red
    exclude_mask: npt.NDArray, optional
        If provided, a bool array of pvals to exclude from plotting
    orange_mask: npt.NDArray, optional
    """
    try:
        pvals = -np.log10(pvals)
    except TypeError:
        # handle Decimals mixed in with np.float64s
        for idx in range(len(pvals)):
            try:
                pvals[idx] = -np.log10(pvals[idx])
            except TypeError:
                pvals[idx] = -np.float64(pvals[idx].log10())
    pvals = pvals.astype(np.float64)
    if exclude_mask is not None:
        positions = positions[exclude_mask]
        pvals = pvals[exclude_mask]
        red_mask = red_mask[exclude_mask]
    if red_mask is not None:
        blue_mask = np.logical_not(red_mask)
        ax.scatter(positions[blue_mask], pvals[blue_mask])
        ax.scatter(positions[red_mask], pvals[red_mask], c="red")
        if orange_mask is not None:
            ax.scatter(positions[orange_mask], pvals[orange_mask], c="orange")
    else:
        ax.scatter(positions, pvals)

def regress_effect(pt: Phenotypes, covars: Genotypes = None):
    """
    Regress out the effect(s) of the covariates on the phenotype

    Parameters
    ----------
    pt: Phenotypes
        A Phenotypes object with only a single phenotype
    covars: Genotypes, optional
        The variables on which to condition the SNPs by encoding them as covariates

    Returns
    -------
    Phenotypes
        The new phenotype values
    """
    resids = Phenotypes(fname=None, log=pt.log)
    resids.samples = pt.samples
    resids.names = pt.names
    if covars is None or not covars.data.shape[1]:
        resids.data = pt.data.copy()
    else:
        resids.data = (
            sm.OLS(pt.data, sm.add_constant(covars.data.sum(axis=2)))
            .fit()
            .resid[:, np.newaxis]
        )
    return resids

def condition_on_variable(gts: npt.NDArray, pt: Phenotypes, covars: Genotypes = None):
    """
    Compute pvals for each SNP or haplotype in a set after conditioning on others

    Parameters
    ----------
    gts: npt.NDArray
        The genotypes of all of the SNPs
    pt: Phenotypes
        A Phenotypes object with only a single phenotype
    covars: Genotypes, optional
        The variables on which to condition the SNPs by encoding them as covariates

    Returns
    -------
    npt.NDArray[float]
        The p-values of all of the input SNPs when conditioned on the covariates
    """
    if covars is None:
        assoc_test = AssocTestSimpleSM()
    else:
        assoc_test = AssocTestSimpleCovariates(covars=covars.data.sum(axis=2))
    return assoc_test.run(gts.sum(axis=2), pt.data[:, 0]).data["pval"]

def condition_on_variable_chunked(
    gts: Genotypes,
    pt: Phenotypes,
    covars: Genotypes = None,
    chunk_size: int = None,
):
    pvals = np.empty(gts.data.shape[1], dtype=object)

    chunks = chunk_size
    if chunks is None or chunks > len(pvals):
        chunks = len(pvals)

    for start in range(0, len(pvals), chunks):
        end = start + chunks
        if end > len(pvals):
            end = len(pvals)
        size = end - start

        pvals[start:end] = condition_on_variable(gts.data[:, start:end], pt, covars)

    return pvals

def extract_snp_matrix(logfile: Path) -> npt.NDArray:

    text = logfile.read_text().splitlines()
    iterations = {}
    current_iter = None
    current_tree = None

    for lineno, line in enumerate(text, start=1):
        ls = line.strip()
        # robustly find "Iteration X: tree Y" even if line is indented
        m = re.search(r"Iteration\s+(\d+)\s*:\s*tree\s+(\d+)", ls, re.I)
        if m:
            current_iter, current_tree = map(int, m.groups())
            continue

        # find label lines (skip lines mentioning 'root')
        if '[label="' in line and 'root' not in line:
            # try to capture up to the literal backslash-n ("\\n") first
            m2 = re.search(r'\[label="([^\\]+?)\\n', line)
            if not m2:
                # fallback: capture up to the next closing quote
                m2 = re.search(r'\[label="([^"]+)"', line)
            if not m2:
                raise ValueError(f"Couldn't parse [label=...] on line {lineno}: {line!r}")
            snp = m2.group(1)

            if current_iter is None or current_tree is None:
                raise ValueError(
                    f"Found SNP label on line {lineno} but no preceding 'Iteration ...: tree ...' line."
                )

            # if there are multiple SNPs in a tree, just choose the first one
            iterations.setdefault(current_iter, {}).setdefault(current_tree, list()).append(snp)

    if not iterations:
        return np.empty((0, 0), dtype=object)

    # remove duplicates in each entry but preserve order of first occurence
    for iteration in iterations:
        for tree in iterations[iteration]:
            iterations[iteration][tree] = tuple(dict.fromkeys(iterations[iteration][tree]))

    # enforce same tree keys for every iteration
    first_it = sorted(iterations.keys())[0]
    expected_keys = sorted(iterations[first_it].keys())
    for it in sorted(iterations.keys()):
        keys = sorted(iterations[it].keys())
        if keys != expected_keys:
            raise ValueError(f"Inconsistent tree keys for iteration {it}: {keys} != {expected_keys}")

    rows = []
    for it in sorted(iterations.keys()):
        # unwrap tuples with a single string but leave multi-element tuples untouched
        row = [
            iterations[it][k][0] if len(iterations[it][k]) == 1 else iterations[it][k]
            for k in expected_keys
        ]
        rows.append(row)

    return np.array(rows, dtype=object)


def flatten_list_of_strings_and_tuples(mixed_list: list):
    flattened_list = []

    for item in mixed_list:
        if isinstance(item, tuple):
            flattened_list.extend(item)
        else:
            flattened_list.append(item)

    return flattened_list


@click.command()
@click.argument("genotypes", type=click.Path(path_type=Path))
@click.argument("phenotype", type=click.Path(path_type=Path))
@click.argument("haplotype", type=click.Path(path_type=Path))
@click.option(
    "-i",
    "--hap-id",
    type=str,
    show_default="all haplotypes",
    help=(
        "A haplotype ID from the .hap file to plot"
        "(ex: '-i H1')."
    ),
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
    show_default="all SNPs",
    help="Only select SNPs with a MAF above this threshold",
)
@click.option(
    "--show-original",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to also depict the original Manhattan plot",
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help=(
        "Perform reading/writing operations in chunks of X variants. "
        "This reduces memory but at the cost of time."
    ),
)
@click.option(
    "--log-file",
    type=Path,
    default=None,
    show_default="no extra plots",
    help=(
        "If provided, we will read trees from this happler log file and attempt to "
        "plot them as additional rows in the conditional regression plot"
    )
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
    genotypes: Path,
    phenotype: Path,
    haplotype: Path,
    hap_id: str = None,
    region: str = None,
    maf: float = None,
    show_original: bool = False,
    chunk_size: int = None,
    log_file: Path = None,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Make Manhattan plots for a set of SNPs for each of three cases:
    1) the haplotype as a covariate
    2) the haplotypes' alleles as (multiple) covariates
    3) each of the haplotypes' alleles as a covariate, where we produce one plot for each allele
    4) the original manhattan plot without any conditioning
    """
    log = getLogger("plot_conditional_regressions", verbosity)

    # load a haplotype and its alleles
    hps = Haplotypes(haplotype, log=log)
    hps.read(haplotypes=(set((hap_id,)) if hap_id is not None else None))
    # get the variants from the haplotype
    # use dict to remove duplicate IDs
    variants = tuple(dict.fromkeys(var.id for hap in hps.data.values() for var in hap.variants))

    # load just the first phenotype
    pts = Phenotypes(phenotype, log=log)
    pts.read()

    # load the SNPs
    gts = GenotypesVCF(genotypes, log=log)
    if genotypes.suffix == ".pgen":
        gts = GenotypesPLINK
    gts = gts(genotypes, log=log, chunk_size=chunk_size)
    gts.read(region=region, samples=set(pts.samples))
    # we need phasing for transform (below)
    gts.check_phase()
    gts.check_missing()
    gts.check_biallelic()
    gts.check_maf(threshold=maf, discard_also=True)
    if not set(hps.data.keys()).isdisjoint(gts.variants["id"]):
        # check that the haplotypes aren't already in the gts and remove them if they are
        gts.subset(variants=tuple(v for v in gts.variants["id"] if v not in hps.data.keys()), inplace=True)
    positions = gts.variants["pos"]
    # reorder the phenotypes and ensure there is only one phenotype
    pts.subset(names=(pts.names[0],), samples=gts.samples, inplace=True)
    pts.standardize()

    # get the haplotype genotypes
    hap_gt = hps.transform(gts)

    # read SNPs from the log file
    other_alleles = np.array([], dtype=object)
    if log_file is None:
        # make the figure
        # set all panels in the same row
        figsize = (FIGSIZE*(len(variants)+2+show_original)/2.5, FIGSIZE)
        fig, axs = plt.subplots(1, 2+len(variants)+show_original, sharey=True, figsize=figsize)
    else:
        show_original = True
        other_alleles = extract_snp_matrix(log_file)

        # make the figure
        # set all panels in the same row
        figsize = (FIGSIZE*(len(variants)+2+show_original)/2.5, FIGSIZE*other_alleles.shape[0])
        fig, axs = plt.subplots(1+other_alleles.shape[0], 2+len(variants)+show_original, sharey=True, figsize=figsize)
        all_axs = axs
        axs = axs[0]

    # Reverse panel placement so plots are filled right-to-left in indexing,
    # which displays left-to-right in the desired reversed order.
    def ax_at(i: int):
        return axs[-(i + 1)]

    # highlight alleles in red
    red_mask = np.zeros(len(gts.variants), dtype=np.bool_)
    for snp in variants:
        red_mask[gts._var_idx[snp]] = True

    log.info("Creating haplotype plot")
    # first, encode the haplotype as covariate
    make_manhattan(
        ax_at(0),
        positions,
        condition_on_variable_chunked(gts, pts, hap_gt, chunk_size=chunk_size),
        red_mask
    )
    ax_at(0).set_title("Haplotype"+("s" if len(hap_gt.variants)-1 else ""))
    log.info("Creating haplotype alleles plot")
    # now, encode the haplotypes' alleles as separate covariates
    covars = gts.subset(variants=variants)
    # exclude the covariates from plotting
    exclude = np.ones(len(gts.variants), dtype=np.bool_)
    for snp in variants:
        exclude[gts._var_idx[snp]] = False
    make_manhattan(
        ax_at(1),
        positions,
        condition_on_variable_chunked(gts, pts, covars, chunk_size=chunk_size),
        red_mask,
        exclude,
    )
    # f-string expressions cannot include a backslash
    f_str_quotation_plural_agh = 's\'' if len(hap_gt.variants)-1 else '\'s'
    ax_at(1).set_title(f"Haplotype{f_str_quotation_plural_agh} Alleles")
    # finally, encode each of the alleles as a covariate in a separate plot
    for idx in range(len(variants)):
        log.info(f"Creating plot #{idx+3}")
        covars = gts.subset(variants=(variants[idx],))
        exclude = np.ones(len(gts.variants), dtype=np.bool_)
        exclude[gts._var_idx[variants[idx]]] = False
        make_manhattan(
            ax_at(idx+2),
            positions,
            condition_on_variable_chunked(gts, pts, covars, chunk_size=chunk_size),
            red_mask,
            exclude,
        )
        ax_at(idx+2).set_title(variants[idx])
    if show_original:
        log.info("Creating original manhattan plot")
        # we're going to have to append the haplotypes to the original SNP genotypes,
        # so we'll have to extend the red mask to include the haplotypes
        # and create an orange mask for the haplotypes themselves
        red_mask = np.append(red_mask, np.zeros(len(hap_gt.variants), dtype=np.bool_))
        orange_mask = np.zeros(red_mask.shape[0], dtype=np.bool_)
        orange_mask[-len(hap_gt.variants):] = 1
        new_gt = Genotypes.merge_variants((gts, hap_gt), fname=None, log=log)
        # now, let's finally plot everything and append the haps to the gts
        make_manhattan(
            ax_at(len(variants)+2),
            new_gt.variants["pos"],
            condition_on_variable_chunked(new_gt, pts, chunk_size=chunk_size),
            red_mask,
            exclude_mask=None,
            orange_mask=orange_mask,
        )
        ax_at(len(variants)+2).set_title("")

    if log_file is not None:
        new_gt.index()
        assert len(hps.data) == 1, "We can only handle one haplotype at a time"
        hap_id = list(hps.data)[0]
        hap_idx = new_gt._var_idx[hap_id]
        flat_other_alleles = other_alleles.flatten()
        for iteration_idx in range(other_alleles.shape[0]):
            covar_variants_for_hap = None
            for tree_idx in range(other_alleles.shape[1]):
                # set up some useful variables for later
                curr_ax = all_axs[iteration_idx+1, tree_idx]
                curr_flat_other_allele_idx = iteration_idx*other_alleles.shape[1] + tree_idx
                start_flat_other_allele_idx = 0
                if iteration_idx:
                    start_flat_other_allele_idx = curr_flat_other_allele_idx - other_alleles.shape[1] + 1
                covar_variants = flat_other_alleles[start_flat_other_allele_idx:curr_flat_other_allele_idx].tolist()
                variant_found = other_alleles[iteration_idx, tree_idx]
                orange_mask = None
                exclude_mask = None
                red_mask = np.zeros(len(new_gt.variants), dtype=np.bool_)
                # did happler find a single SNP or a haplotype in this tree?
                if isinstance(variant_found, tuple):
                    covar_variants_for_hap = covar_variants
                    # verify that the haplotype matches the one that was provided
                    assert sorted(variant_found) == sorted(variants), "Encountered a haplotype which was not in the .hap file"
                    # highlight the haplotype that we found in the tree as orange
                    orange_mask = np.zeros(len(new_gt.variants), dtype=np.bool_)
                    orange_mask[hap_idx] = True
                    # now, also highlight the individual alleles
                    for hap_snp in variant_found:
                        red_mask[new_gt._var_idx[hap_snp]] = True
                    variant_found = hap_id
                else:
                    # exclude the haplotype from plotting
                    exclude_mask = np.ones(len(new_gt.variants), dtype=np.bool_)
                    exclude_mask[hap_idx] = False
                    # highlight the variant that we found in the tree as red
                    red_mask[new_gt._var_idx[variant_found]] = True
                # is one of the covar_variants a haplotype?
                has_hap = [isinstance(item, tuple) for item in covar_variants]
                if sum(has_hap):
                    try:
                        assert variants in covar_variants
                    except:
                        assert variants[::-1] in covar_variants
                        covar_variants[covar_variants.index(variants[::-1])] = hap_id
                    else:
                        covar_variants[covar_variants.index(variants)] = hap_id
                # now we can finally regress things out
                if covar_variants:
                    # set title to indicate the variants that we are conditioning on
                    curr_ax.set_title("\n".join(covar_variants))
                    covars = new_gt.subset(variants=set(covar_variants))
                    resids = regress_effect(pts, covars)
                else:
                    resids = pts
                curr_ax.text(.5,.96, variant_found, fontsize=9, horizontalalignment='center', transform=curr_ax.transAxes)
                make_manhattan(
                    curr_ax,
                    new_gt.variants["pos"],
                    condition_on_variable_chunked(new_gt, resids, chunk_size=chunk_size),
                    exclude_mask=exclude_mask,
                    red_mask=red_mask,
                    orange_mask=orange_mask,
                )
            for axes_idx in range(all_axs.shape[1]-1, other_alleles.shape[1]-1, -1):
                curr_ax = all_axs[iteration_idx+1, axes_idx]
                if covar_variants_for_hap is None or axes_idx != (all_axs.shape[1] - 1):
                    curr_ax.remove()
                    continue
                # If this is the last axis and we saw a haplotype, let's make a midway plot
                assert len(variants) == 2, "This will only work for a haplotype with two alleles"
                curr_ax.set_title(f"midway {hap_id}:\n"+"\n".join(covar_variants_for_hap))
                covars = gts.subset(variants=set(covar_variants_for_hap))
                resids = regress_effect(pts, covars)
                # now, transform all of the SNPs by the haplotype up until the target_variant
                gts_tsfm = Genotypes(fname=None, log=log)
                haptools_haplotype = copy.deepcopy(list(hps.data.values())[0])
                target_variant = haptools_haplotype.variants[-1]
                target_allele = gts.variants[gts._var_idx[target_variant.id]]["alleles"][0]
                target_allele = int(target_variant.allele != target_allele)
                haptools_haplotype.variants = haptools_haplotype.variants[:-1]
                hp = HapplerHaplotype.from_haptools_haplotype(haptools_haplotype, gts)
                gts_tsfm.variants = np.delete(gts.variants, hp.node_indices)
                gts_tsfm.samples = gts.samples
                gts_tsfm.data = hp.transform(gts, target_allele)
                gts_tsfm.index()
                # also make a red mask for the variant
                red_mask = np.zeros(len(gts_tsfm.variants), dtype=np.bool_)
                red_mask[gts_tsfm._var_idx[target_variant.id]] = True
                make_manhattan(
                    curr_ax,
                    gts_tsfm.variants["pos"],
                    condition_on_variable_chunked(gts_tsfm, resids, chunk_size=chunk_size),
                    exclude_mask=None,
                    red_mask=red_mask,
                    orange_mask=None,
                )

    # now, tidy up and save the plot
    log.info("Writing out plot")
    fig.supylabel("-log10(pval)")
    fig.supxlabel("Chromosomal Position")
    fig.suptitle("Manhattan plots when conditioning on...")
    fig.tight_layout()
    fig.savefig(output)

if __name__ == "__main__":
    main()
