#!/usr/bin/env python
import csv
import copy
import pickle
from pathlib import Path
from logging import Logger
from dataclasses import dataclass, field

import click
import matplotlib
import numpy as np
matplotlib.use('Agg')
import numpy.typing as npt
import statsmodels.api as sm
import matplotlib.pyplot as plt

from haptools.logging import getLogger
from haptools.ld import pearson_corr_ld
from haptools.data import Phenotypes, GenotypesVCF, GenotypesPLINK, Haplotypes
from haptools.data import Extra, Haplotype as HaplotypeBase, Variant as VariantBase

from snakemake_io import glob_wildcards


@dataclass
class HapplerVariant(VariantBase):
    """
    A variant allele with sufficient fields for happler
    Properties and functions are shared with the base Variant object, "VariantBase"
    """

    score: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(Extra("score", ".2f", "BIC assigned to this variant"),),
    )


@dataclass
class HapplerHaplotype(HaplotypeBase):
    """
    A haplotype with sufficient fields for happler
    Properties and functions are shared with the base Haplotype object, "HaplotypeBase"
    """

    beta: float
    pval: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(
            Extra("beta", ".2f", "Effect size in linear model"),
            Extra("pval", ".2f", "-log(pval) in linear model"),
        ),
    )

def standardize(data):
    """
    Standardize a matrix so it has a mean of 0 and a stdev of 1
    """
    if len(data.shape) <= 1:
        raise ValueError("The data property must have a 2D shape.")
    std = np.std(data, axis=0)
    new_data = (data - np.mean(data, axis=0)) / std
    # for phenotypes where the stdev is 0, just set all values to 0 instead of nan
    zero_elements = std == 0
    new_data[:, zero_elements] = np.zeros(
        (new_data.shape[0], np.sum(zero_elements))
    )
    return new_data


def get_beta_and_rsquared(result: sm.regression.linear_model.RegressionResults):
    return result.params[0], result.rsquared


def compute_explained_variance(gt: npt.NDArray, pt: npt.NDArray):
    """
    Compute explained variance for each SNP or haplotype in a set
    Also compute R-squared values

    The explained variance is beta^2 in a linear model y = bx + e when x and y have
    been standardized to mean 0 and stdev 1

    Parameters
    ----------
    gt: npt.NDArray
        The genotypes of the SNP (or haplotype) as a 2D array of shape:
        num_samples by num_variants
    pt: npt.NDArray
        The phenotype as a 1D array

    Returns
    -------
    float
        The explained variance of the SNPs/haplotypes and the R squared values
    """
    # standardize the phenotypes and genotypes
    pt = standardize(pt[:, np.newaxis]).flatten()
    gt = standardize(gt)
    # compute the betas of each SNP by fitting a linear model y = beta * x + e
    betas_rsquared = np.array([
        get_beta_and_rsquared(sm.OLS(pt, gt[:, snp_idx].flatten()).fit())
        for snp_idx in range(gt.shape[1])
    ])
    # compute explained variance by squaring the betas
    betas_rsquared[:, 0] = betas_rsquared[:, 0]**2
    return betas_rsquared


def compute_multisnp_rsquared(gt: npt.NDArray, pt: npt.NDArray):
    """
    Compute the R-squared value for SNPs (from a single haplotype) in a multiple linear
    regression

    Parameters
    ----------
    gt: npt.NDArray
        The genotypes of the SNPs in a haplotype as a 2D array of shape:
        num_samples by num_variants
    pt: npt.NDArray
        The phenotype as a 1D array

    Returns
    -------
    float
        The R-squared value of the SNPs from a haplotype
    """
    # standardize the phenotypes and genotypes
    pt = standardize(pt[:, np.newaxis]).flatten()
    gt = standardize(gt)
    # fitting a linear model y = beta1 * x1 + beta2 * x2 + ... + e
    return sm.OLS(pt, gt).fit().rsquared


def compute_summary_stats(gt: npt.NDArray, pt: npt.NDArray):
    """
    Compute summary stats for SNPs (from a single haplotype) in a multiple linear
    regression

    Parameters
    ----------
    gt: npt.NDArray
        The genotypes of a variant in a 1D array of length num_samples
    pt: npt.NDArray
        The phenotype as a 1D array

    Returns
    -------
    tuple[float, float, float]
        1. Effect size
        2. P-value
        3. BIC
    """
    # standardize the phenotypes and genotypes
    pt = standardize(pt[:, np.newaxis]).flatten()
    gt = standardize(gt)
    # fitting a linear model y = beta1 * x1
    fit = sm.OLS(pt, gt).fit()
    return (fit.params[0], fit.pvalues[0], fit.bic)


def load_data(
    gts: Path,
    hps: Path,
    pts: Path,
    log: Logger = None
):
    """
    Load gts, hps, and pts properly

    Parameters
    ----------
    gts: Path
        The path to a PGEN file containing genotypes for all haplotypes and their SNPs
    hps: Path
        The path to a .hap file containing a set of haplotypes
    pts: Path
        The path to a pheno file containing the phenotypes
    log: Logger, optional
        A logging object to write any debugging and error messages
    
    Returns
    -------
    tuple[gts, hps, pts]
    """
    # load the phenotypes
    pts = Phenotypes(pts, log=log)
    pts.read()

    # load the haplotypes
    hps_path = hps
    hps = Haplotypes(hps, haplotype=HapplerHaplotype, variant=HapplerVariant, log=log)
    hps.read()
    if not len(hps.data):
        return dict()

    # which variants do we need?
    snp_variants = {v.id for hap in hps.data.values() for v in hap.variants}
    variants = snp_variants | hps.data.keys()

    # what region should we use?
    # we can figure it out by looking at the haps in the .hap file
    chrom = next(iter(hps.data.values())).chrom
    min_pos = min(v.start for hap in hps.data.values() for v in hap.variants)
    max_pos = max(v.end for hap in hps.data.values() for v in hap.variants)
    region = chrom + ":" + str(min_pos) + "-" + str(max_pos)

    # load the SNP and hap genotypes
    gts_path = gts
    gts = GenotypesPLINK(fname=gts, log=log)
    gts.read(variants=variants, region=region, samples=set(pts.samples))
    gts.check_phase()
    gts.check_missing()
    gts.check_biallelic()
    gts.index()
    assert snp_variants.issubset(gts.variants["id"]), f"Couldn't find all variants for {hps_path} in {gts_path}"
    # check that the hps IDs are in there too
    if not all(h in gts.variants["id"] for h in hps.data.keys()):
        gts = GenotypesPLINK.merge_variants((gts, hps.transform(gts)), fname=gts.fname)
    assert len(gts.variants) == len(variants)
    pts.subset(samples=gts.samples, inplace=True)

    return gts, hps, pts

def get_explained_variances(
    gts: Path,
    hps: Path,
    pts: Path,
    log: Logger = None,
):
    """
    Compute explained variance for the haplotypes in a .hap file and each haplotypes'
    SNPs

    Parameters
    ----------
    gts: Path
        The path to a PGEN file containing genotypes for all haplotypes and their SNPs
    hps: Path
        The path to a .hap file containing a set of haplotypes
    pts: Path
        The path to a pheno file containing the phenotypes
    log: Logger, optional
        A logging object to write any debugging and error messages

    Returns
    -------
    dict[str, tuple[float, float, float, float]]
        The dict is keyed by each haplotype's ID and has the following values:
        1. explained variance for the haplotype
        2. R-squared for the haplotype
        3. explained variance for the haplotype's SNPs
        4. R-squared for the haplotype's SNPs
    """
    gts, hps, pts = load_data(gts, hps, pts)

    # compute explained variance for each SNP and haplotype
    raw_explained_variances = dict(zip(gts.variants["id"], compute_explained_variance(
        gts.data.sum(axis=2), pts.data[:, 0],
    )))

    # create a dictionary mapping hap IDs to two-element tuples containing
    # 1) explained variance for the haplotype
    # 2) R-squared for the haplotype
    # 3) explained variance for the haplotype's SNPs
    # 4) R-squared for the haplotype's SNPs
    vals = {}
    for hp in hps.data.values():
        sum_of_SNPs = sum(raw_explained_variances[v.id][0] for v in hp.variants)
        multi_rsquared = compute_multisnp_rsquared(
            gts.subset(variants=tuple(v.id for v in hp.variants)).data.sum(axis=2),
            pts.data[:, 0],
        )
        vals[hp.id] = (
            *raw_explained_variances[hp.id],
            sum_of_SNPs,
            multi_rsquared,
        )

    return vals


def get_metrics(
    gts: Path,
    hps: Path,
    pts: Path,
    log: Logger = None
):
    """
    Get metrics for the haplotypes and variants in this .hap file

    Parameters
    ----------
    gts: Path
        The path to a PGEN file containing genotypes for all haplotypes and their SNPs
    hps: Path
        The path to a .hap file containing a set of haplotypes
    pts: Path
        The path to a pheno file containing the phenotypes
    log: Logger, optional
        A logging object to write any debugging and error messages

    Returns
    -------
    dict[str, tuple[float, float, float, float]]
        The dict is keyed by each haplotype's ID and has the following values:
        0. LD between the alleles of the haplotype
        1. allele frequency for the haplotype
        2. effect size for the haplotype
        3. p-value for the haplotype
        4. BIC value for the haplotype
        5. BIC value for all alleles of the haplotype (independently)
        6. the number of alleles in the haplotype
        7. allele frequency for the haplotype's SNPs
        8. effect size for the haplotype's SNPs
        9. p-value for the haplotype's SNPs
        10. BIC value for the haplotype's SNPs
        11. delta BIC values for the haplotype's SNPs
        12. Original BIC value for the haplotype's SNPs (from the .hap file)
        13. Original delta BIC values for the haplotype's SNPs (from the .hap file)
    """
    gts, hps, pts = load_data(gts, hps, pts)

    afs = dict(zip(gts.variants['id'], gts.check_maf()))

    summary_stats = {}
    for variant in gts.variants["id"]:
        summary_stats[variant] = compute_summary_stats(
            gts.subset(variants=(variant,)).data.sum(axis=2),
            pts.data[:, 0],
        )

    vals = {}
    for hp in hps.data.values():
        hp_vars = tuple(variant.id for variant in hp.variants)
        hp_vars_alls = {variant.id: variant.allele for variant in hps.data[hp.id].variants}
        hp_vars_gt = gts.subset(variants=hp_vars)
        # do we need to flip the allele?
        hp_vars_alls = {
            als[0]: als[1][0] == hp_vars_alls[als[0]]
            for als in hp_vars_gt.variants[["id","alleles"]]
        }
        if len(hp_vars) > 2:
            hp_vars_ld = 0
        elif len(hp_vars) == 2:
            hp_vars_ld = pearson_corr_ld(
                (hp_vars_gt.data[:,0] == list(hp_vars_alls.values())[0]).sum(axis=1),
                (hp_vars_gt.data[:,1] == list(hp_vars_alls.values())[1]).sum(axis=1)
            )
        else:
            hp_vars_ld = float("inf")
        indep_bic = compute_summary_stats(
            gts.subset(variants=hp_vars).data.sum(axis=2),
            pts.data[:, 0],
        )[2]
        hp_al_afs = np.empty(len(hp_vars))
        hp_al_effects = np.empty(len(hp_vars))
        hp_al_pvals = np.empty(len(hp_vars))
        hp_al_bics = np.empty(len(hp_vars))
        hp_al_deltas = np.empty(len(hp_vars))
        for hp_al_idx in range(len(hp_vars)):
            hp_cp = copy.deepcopy(hp)
            hp_cp.variants = hp.variants[:hp_al_idx+1]
            hp_cp_gt = hp_cp.transform(gts).sum(axis=1)[:, np.newaxis]
            hp_al_summary_stats = compute_summary_stats(hp_cp_gt, pts.data[:, 0])

            hp_al_afs[hp_al_idx] = hp_cp_gt.sum() / (2*len(hp_cp_gt))
            hp_al_effects[hp_al_idx] = hp_al_summary_stats[0]
            hp_al_pvals[hp_al_idx] = hp_al_summary_stats[1]
            hp_al_bics[hp_al_idx] = hp_al_summary_stats[2]
        hp_al_bics_og = np.array([v.score for v in hps.data[hp.id].variants])
        hp_al_deltas = hp_al_bics[:-1] - hp_al_bics[1:]
        hp_al_deltas_og = hp_al_bics_og[:-1] - hp_al_bics_og[1:]
        vals[hp.id] = (
            hp_vars_ld,
            afs[hp.id],
            summary_stats[hp.id][0],
            summary_stats[hp.id][1],
            summary_stats[hp.id][2],
            indep_bic,
            len(hp_vars),
            list(hp_al_afs),
            list(hp_al_effects),
            list(hp_al_pvals),
            list(hp_al_bics),
            list(hp_al_deltas),
            list(hp_al_bics_og),
            list(hp_al_deltas_og),
        )

    return vals


@click.command()
@click.argument("genotypes", type=click.Path(path_type=Path))
@click.argument("phenotypes", type=click.Path(path_type=Path))
@click.argument("haplotypes", type=click.Path(path_type=Path))
@click.option(
    "-s",
    "--subset",
    type=click.Path(path_type=Path),
    default=None,
    show_default="no subsetting",
    help="A .txt file containing a subset of haplotype files to consider",
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
    phenotypes: Path,
    haplotypes: Path,
    subset: Path = None,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Plot explained variation of haplotypes vs SNPs for a bunch of loci

    Each locus is inferred from brace expressions injected into the paths to the files
    For example, a path like "out/{region}/happler/hap/happler.hap"
    will infer that {region} is a locus. Wildcards will be expanded across all inputs.
    """
    log = getLogger("plot_variance_explained", verbosity)

    # extract parameters and parameter values by globbing wildcards
    # params will be a dictionary mapping parameter names to lists of values
    params = dict(glob_wildcards(haplotypes)._asdict())
    dtypes = {k: "U30" for k in params.keys()}
    # convert the dictionary to a numpy mixed dtype array
    params = np.array(list(zip(*params.values())), dtype=list(dtypes.items()))

    get_hap_fname = lambda path, param_set: Path(str(path).format(**dict(zip(dtypes.keys(), param_set))))

    # figure out which params to use based on the subset file
    if subset is not None:
        with open(subset, "r") as f:
            subset = set(map(Path, f.read().splitlines()))
        subset_mask = [
            get_hap_fname(haplotypes, params[idx]) in subset
            for idx in range(len(params))
        ]
        params = params[subset_mask]
        if not len(params):
            log.warning("All files were removed after subsetting")

    # compute explained variance for each haplotype and its SNPs
    vals = [
        (params[idx], hp_id, hap)
        for idx in range(len(params))
        for hp_id, hap in get_explained_variances(
            get_hap_fname(genotypes, params[idx]),
            get_hap_fname(haplotypes, params[idx]),
            get_hap_fname(phenotypes, params[idx]),
            log=log
        ).items()
    ]

    # compute miscellaneous metrics
    other_vals = [
        (params[idx], hp_id, hap)
        for idx in range(len(params))
        for hp_id, hap in get_metrics(
            get_hap_fname(genotypes, params[idx]),
            get_hap_fname(haplotypes, params[idx]),
            get_hap_fname(phenotypes, params[idx]),
            log=log
        ).items()
    ]

    # this 2D array should have two * 2 columns: 1) the haplotype and 2) its SNPs
    # and should have as many rows as there are haplotypes among all of the loci
    # (note that some loci may have multiple haplotypes so we adjust 'params' accordingly)
    params, hp_ids, vals = np.array([v[0] for v in vals]), [v[1] for v in vals], np.array([v[2] for v in vals])
    explained_variances = vals[:, (0, 2)]
    rsquareds = vals[:, (1, 3)]

    if np.any(explained_variances[:, 0] > 1):
        log.error(
            "Some of the explained variances are greater than 1! Check that nothing "
            "went wrong."
        )

    params1, hp_ids1, other_vals = np.array([v[0] for v in other_vals]), [v[1] for v in other_vals], np.array([v[2] for v in other_vals], dtype=object)
    assert (params == params1).all()
    assert hp_ids1 == hp_ids
    lds = other_vals[:, 0]
    hap_afs = other_vals[:, 1]
    hap_betas = other_vals[:, 2]
    hap_pvals = other_vals[:, 3]
    hap_bics = other_vals[:, 4]
    indep_bics = other_vals[:, 5]
    num_alleles = other_vals[:, 6]
    allele_afs = other_vals[:, 7]
    allele_effects = other_vals[:, 8]
    allele_pvals = other_vals[:, 9]
    allele_bics = other_vals[:, 10]
    allele_deltas = other_vals[:, 11]
    allele_bics_og = other_vals[:, 12]
    allele_deltas_og = other_vals[:, 13]

    # how good are we doing?
    percent_success = 100*(
        explained_variances[:, 1] < explained_variances[:, 0]
    ).sum()/explained_variances.shape[0]
    log.info(
        f"{percent_success:.3f}% of haplotypes explain more phenotypic variation than"
        " their SNPs"
    )
    percent_success = 100*(
        explained_variances[:, 1] == explained_variances[:, 0]
    ).sum()/explained_variances.shape[0]
    log.info(
        f"{percent_success:.3f}% of haplotypes explain as much phenotypic variation as"
        " their SNPs"
    )
    # how good are we doing?
    percent_success = 100*(
        rsquareds[:, 1] < rsquareds[:, 0]
    ).sum()/rsquareds.shape[0]
    log.info(
        f"{percent_success:.3f}% of haplotypes explain more R-squared than"
        " their SNPs"
    )
    percent_success = 100*(
        rsquareds[:, 1] == rsquareds[:, 0]
    ).sum()/rsquareds.shape[0]
    log.info(
        f"{percent_success:.3f}% of haplotypes explain as much R-squared as"
        " their SNPs"
    )

    max_ev_val = explained_variances.max()
    max_r2_val = rsquareds.max()

    f, (ax1, ax2) = plt.subplots(1, 2)

    with open(output.with_suffix(".pickle"), "wb") as picklef:
        pickle.dump((
            params,
            explained_variances,
            rsquareds,
            lds,
            hap_afs,
            hap_betas,
            hap_pvals,
            hap_bics,
            indep_bics,
            num_alleles,
            allele_afs,
            allele_effects,
            allele_pvals,
            allele_bics,
            allele_deltas,
            allele_bics_og,
            allele_deltas_og,
        ), picklef)

    with open(output.with_suffix(".tsv"), 'w', newline='') as tsvfile:
        tsv_writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        tsv_writer.writerow([
            "locus",
            "hap_exp_var",
            "alleles_exp_var",
            "hap_r2",
            "alleles_r2",
            "hap_over_alleles_r2",
            "allele_ld",
            "hap_afs",
            "hap_betas",
            "hap_pvals",
            "hap_bics",
            "indep_bics",
            "num_alleles",
            "allele_afs",
            "allele_effects",
            "allele_pvals",
            "allele_bics",
            "allele_deltas",
            "allele_bics_og",
            "allele_deltas_og",
        ])
        if ("locus" in params.dtype.names) and ("gene" in params.dtype.names):
            name = lambda i: i["locus"]+":"+i["gene"]
        else:
            name = lambda i: ":".join(i)
        for i in range(explained_variances.shape[0]):
            tsv_writer.writerow([
                name(params[i])+":"+hp_ids[i],
                explained_variances[i, 0],
                explained_variances[i, 1],
                rsquareds[i, 0],
                rsquareds[i, 1],
                rsquareds[i,0]/rsquareds[i,1],
                lds[i],
                hap_afs[i],
                hap_betas[i],
                hap_pvals[i],
                hap_bics[i],
                indep_bics[i],
                num_alleles[i],
                ",".join(map(str,allele_afs[i])),
                ",".join(map(str,allele_effects[i])),
                ",".join(map(str,allele_pvals[i])),
                ",".join(map(str,allele_bics[i])),
                ",".join(map(str,allele_deltas[i])),
                ",".join(map(str,allele_bics_og[i])),
                ",".join(map(str,allele_deltas_og[i])),
            ])

    ax1.scatter(explained_variances[:, 0], explained_variances[:, 0]/explained_variances[:, 1])
    ax1.axline([0, 1], [max_ev_val, 1])
    ax1.set_title("Variance Explained")
    ax1.set_xlabel("Haplotype")
    ax1.set_ylabel("Haplotype / Haplotype's SNPs")

    ax2.scatter(rsquareds[:, 0], rsquareds[:, 0]/rsquareds[:, 1])
    ax2.axline([0, 1], [max_r2_val, 1])
    ax2.set_title("R-Squared")
    ax2.set_xlabel("Haplotype")
    ax2.set_ylabel("Haplotype / Haplotype's SNPs")

    f.set_size_inches(10, 5)
    plt.tight_layout()
    plt.savefig(output)


if __name__ == "__main__":
    main()
