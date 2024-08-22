#!/usr/bin/env python
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
from haptools.data import Phenotypes, GenotypesVCF, GenotypesPLINK, Haplotypes

from snakemake_io import glob_wildcards


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
    # TODO: should we return adjusted r-squared instead?
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


def get_explained_variances(
    gts: Path,
    hps: Path,
    pts: Path,
    log: Logger = None
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
    dict[str, tuple[float, float]]
        Explained variances of 1) each haplotype in the .hap file and 2) the
        haplotype's SNPs. The dict is keyed by each haplotype's ID
    """
    # load the phenotypes
    pts = Phenotypes(pts, log=log)
    pts.read()

    # load the haplotypes
    hps = Haplotypes(hps, log=log)
    hps.read()

    # which variants do we need?
    variants = {v.id for hap in hps.data.values() for v in hap.variants}
    variants.update(hps.data.keys())

    # what region should we use?
    # we can figure it out by looking at the haps in the .hap file
    chrom = next(iter(hps.data.values())).chrom
    min_pos = min(v.start for hap in hps.data.values() for v in hap.variants)
    max_pos = max(v.end for hap in hps.data.values() for v in hap.variants)
    region = chrom + ":" + str(min_pos) + "-" + str(max_pos)

    # load the SNP and hap genotypes
    gts = GenotypesPLINK(fname=gts, log=log)
    gts.read(variants=variants, region=region, samples=set(pts.samples))
    gts.check_phase()
    gts.check_missing()
    gts.check_biallelic()
    gts.index()
    pts.subset(samples=gts.samples, inplace=True)

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

    # compute explained variance for each haplotype and its SNPs
    # this 2D array should have two*2 columns: 1) the haplotype and 2) its SNPs
    # and should have as many rows as there are haplotypes among all of the loci
    # (note that some loci may have multiple haplotypes)
    vals = np.array([
        hap
        for idx in range(len(params))
        for hap in get_explained_variances(
            get_hap_fname(genotypes, params[idx]),
            get_hap_fname(haplotypes, params[idx]),
            get_hap_fname(phenotypes, params[idx]),
            log=log
        ).values()
    ])
    explained_variances = vals[:, (0, 2)]
    rsquareds = vals[:, (1, 3)]

    if np.any(explained_variances > 1):
        log.error(
            "Some of the explained variances are greater than 1! Check that nothing "
            "went wrong."
        )

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

    ax1.scatter(explained_variances[:, 1], explained_variances[:, 0])
    ax1.axline([0, 0], [max_ev_val, max_ev_val])
    ax1.set_title("Variance Explained")
    ax1.set_xlabel("Haplotype's SNPs")
    ax1.set_ylabel("Haplotype")

    ax2.scatter(rsquareds[:, 1], rsquareds[:, 0])
    ax2.axline([0, 0], [max_r2_val, max_r2_val])
    ax2.set_title("R-Squared")
    ax2.set_xlabel("Haplotype's SNPs")
    ax2.set_ylabel("Haplotype")

    f.set_size_inches(10, 5)
    plt.tight_layout()
    plt.savefig(output)


if __name__ == "__main__":
    main()
