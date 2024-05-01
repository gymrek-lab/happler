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
    output: Path = Path("/dev/stdout"),
    verbosity: str = "DEBUG",
):
    """
    Regress out covariates (in .covar format) from phenotypes (in .pheno format)
    """
    log = getLogger("compute_residuals", verbosity)

    log.info("Loading data")
    cv = Covariates(covariates, log=log)
    pt = Phenotypes(phenotypes, log=log)
    pt.read()
    cv.read(samples=set(pt.samples))
    cv.subset(samples=pt.samples, inplace=True)
    assert pt.samples == cv.samples

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
