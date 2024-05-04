#!/usr/bin/env python
from pathlib import Path
from logging import Logger

import click
import numpy as np
import numpy.typing as npt
import statsmodels.api as sm

from haptools.logging import getLogger
from haptools.data import Phenotypes, Covariates


def compute_residuals(pt: npt.NDArray, cv: npt.NDArray):
    """
    Regress covariates out of phenotypes

    Parameters
    ----------
    pt: npt.NDArray
        The phenotypes as a 2D array (shape: samples x phenotypes)
    cv: npt.NDarray
        The covariates as a 2D array (shape: samples x covariates)

    Returns
    -------
    npt.NDArray
        A 2D array (shape: samples x residuals) containing the residuals of a model
        y ~ X where y refers to the phenotypes and X refers to the covariates
    """
    # compute the residuals of each phenotype by fitting a linear model y = X + e
    resids = np.array([
        sm.OLS(pt[:, idx], cv).fit().resid
        for idx in range(pt.shape[1])
    ]).T
    return resids


@click.command()
@click.argument("phenotypes", type=click.Path(path_type=Path))
@click.argument("covariates", type=click.Path(path_type=Path))
@click.option(
    "-e",
    "--extra-covariates",
    type=click.Path(path_type=Path),
    default=None,
    show_default=True,
    help="The path to an extra .covar file to use",
)
@click.option(
    "-s",
    "--samples-to-remove",
    type=click.Path(path_type=Path),
    default=None,
    show_default=True,
    help="The path to a file containing a list of sample IDs to remove",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="The path to a .pheno file into which to write the residuals",
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
    phenotypes: Path,
    covariates: Path,
    extra_covariates: Path = None,
    samples_to_remove: Path = None,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Regress out covariates (in .covar format) from phenotypes (in .pheno format)

    Covariates are expected to have at least as many samples as the phenotypes
    But --extra-covariates may have fewer or more samples, in which case we will output
    the intersection of the set of samples
    """
    log = getLogger("compute_residuals", verbosity)

    log.info("Loading data")
    pt = Phenotypes(phenotypes, log=log)
    pt.read()
    if samples_to_remove is not None:
        # figure out which samples to remove
        with open(samples_to_remove, "r") as myfile:
            samples_to_remove = set(myfile.read().splitlines())
        samples_to_keep = tuple(
            samp for samp in pt.samples
            if samp not in samples_to_remove
        )
        num_removed = len(pt.samples) - len(samples_to_keep)
        log.info(
            f"Removing {num_removed} of {len(samples_to_remove)} withdrawn samples and"
            f" keeping {len(samples_to_keep)}"
        )
        pt.subset(samples=samples_to_keep, inplace=True)
    samples = set(pt.samples)
    cv = Covariates(covariates, log=log)
    cv.read(samples=samples)
    cv.subset(samples=pt.samples, inplace=True)
    assert pt.samples == cv.samples

    if extra_covariates is not None:
        log.info("Reading extra covariates")
        ecv = Covariates(extra_covariates, log=log)
        ecv.read(samples=samples)
        if len(pt.samples) > len(ecv.samples):
            samples_to_keep = set(ecv.samples)
            samples_to_keep = tuple(
                samp for samp in pt.samples if samp in samples_to_keep
            )
            pt.subset(samples=samples_to_keep, inplace=True)
            cv.subset(samples=pt.samples, inplace=True)
            log.info(
                f"Dropped {len(samples)-len(pt.samples)} of {len(samples)} samples "
                f"that were absent in --extra-covariates. Retaining {len(pt.samples)}"
            )
            samples = set(pt.samples)
        ecv.subset(samples=pt.samples, inplace=True)
        cv.names = cv.names + ecv.names
        log.debug("Appending extra covariates")
        cv.data = np.concatenate((cv.data, ecv.data), axis=1)
        cv._samp_idx = None
        cv._name_idx = None

    log.info("Standardizing input")
    pt.standardize()
    cv.standardize()

    log.info("Initializing output writer")
    resids = Phenotypes(output, log=log)
    resids.samples = pt.samples
    resids.names = pt.names

    log.info("Computing residuals")
    resids.data = compute_residuals(pt.data, cv.data)

    log.info("Standardizing and writing output")
    resids.standardize()
    resids.write()


if __name__ == "__main__":
    main()
