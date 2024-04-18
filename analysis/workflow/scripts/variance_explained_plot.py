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


def compute_explained_variance(beta: float, gt: npt.NDArray, pt: npt.NDArray, override: bool = True):
    """
    compute explained variance for a SNP or haplotype

    If beta is a 1D array instead of a single value, gt should be a 2D array of
    num_samples by num_variants

    Parameters
    ----------
    beta: float
        The effect size of the SNP or haplotype
    gt: npt.NDArray
        The genotypes of the SNP (or haplotype) as a 1D array
    pt: npt.NDArray
        The phenotype as a 1D array
    override: bool, optional
        If True, recompute beta and override it
    
    Returns
    -------
    float
        The explained variance of the SNP or haplotype
    """
    if override:
        betas = np.empty(len(beta), dtype=np.float64)
        pt = standardize(pt[:, np.newaxis]).flatten()
        gt = standardize(gt)
        for j in range(len(beta)):
            model = sm.OLS(pt, gt[:, j].flatten())
            betas[j] = (model.fit().params[0])**2
        return betas
    else:
        return (beta**2) * np.var(gt, axis=0) / np.var(pt)


def get_explained_variances(
    gts: Path,
    hps: Path,
    pts: Path,
    linear: Path,
    region: str = None,
    log: Logger = None
):
    """return explained variance for a hap file and the alleles in it"""
    # fix the region
    region = region.replace("_", ":")

    # load the haplotypes
    hps = Haplotypes(hps, log=log)
    hps.read()

    # which variants do we need?
    variants = {v.id for hap in hps.data.values() for v in hap.variants}
    variants.update(hps.data.keys())

    # load the SNP and hap genotypes
    gts = GenotypesPLINK(fname=gts, log=log)
    gts.read(variants=variants, region=region)
    gts.check_phase()
    gts.check_missing()
    gts.check_biallelic()
    gts.index()

    # load the phenotypes
    pts = Phenotypes(pts, log=log)
    pts.read()

    # load the linear files
    ids = [None]*len(gts.variants)
    betas = np.empty(len(gts.variants), dtype=np.float64)
    with open(linear, 'r') as linear_file:
        # Read the first line to find column indices
        header = linear_file.readline().split("\t")
        id_index = header.index('ID')
        beta_index = header.index('BETA')
        # Read the rest of the lines
        idx = 0
        for line in linear_file:
            parts = line.split("\t")
            if parts[id_index] in variants:
                ids[idx] = parts[id_index]
                betas[idx] = parts[beta_index]
                idx += 1

    # reorder betas in case their ordering is different than in gts.variants
    betas = betas[[gts._var_idx[v] for v in ids]]

    # compute explained variance for each SNP and haplotype
    raw_explained_variances = dict(zip(gts.variants["id"], compute_explained_variance(
        betas, gts.data.sum(axis=2), pts.data[:, 0],
    )))

    # create a dictionary mapping hap IDs to two-element tuples containing
    # explained variances for 1) the haplotype and 2) the sum of the haplotype's SNPs
    explained_variances = {}
    for hp in hps.data.values():
        sum_of_SNPs = sum(raw_explained_variances[v.id] for v in hp.variants)
        if sum_of_SNPs > 1:
            breakpoint()
        explained_variances[hp.id] = (raw_explained_variances[hp.id], sum_of_SNPs)

    return explained_variances.values()


@click.command()
@click.argument("genotypes", type=click.Path(path_type=Path))
@click.argument("phenotypes", type=click.Path(path_type=Path))
@click.argument("haplotypes", type=click.Path(path_type=Path))
@click.argument("linears", type=click.Path(path_type=Path))
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
    linears: Path,
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
    # this 2D array should have two columns: 1) the haplotype and 2) its SNPs
    # and should have as many rows as there are haplotypes among all of the loci
    # (note that some loci may have multiple haplotypes)
    explained_variances = np.array([
        hap
        for idx in range(len(params))
        for hap in get_explained_variances(
            get_hap_fname(genotypes, params[idx]),
            get_hap_fname(haplotypes, params[idx]),
            get_hap_fname(phenotypes, params[idx]),
            get_hap_fname(linears, params[idx]),
            region=(params[idx]["locus"] if "locus" in dtypes else None),
            log=log
        )
    ])

    percent_success = (
        explained_variances[:, 1] < explained_variances[:, 0]
    ).sum()/explained_variances.shape[0]
    print(percent_success)

    max_val = explained_variances.max()

    plt.scatter(explained_variances[:, 1], explained_variances[:, 0])
    plt.axline([0, 0], [max_val, max_val])
    plt.title("Variance Explained")
    plt.xlabel("Haplotype's SNPs")
    plt.ylabel("Haplotype")
    plt.savefig(output)


if __name__ == "__main__":
    main()
