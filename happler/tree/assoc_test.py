from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass
from abc import ABC, abstractmethod

import numpy as np
from scipy import stats
import numpy.typing as npt
import statsmodels.api as sm


# We declare this class to be a dataclass to automatically define __init__ and a few
# other methods
@dataclass
class AssocResults:
    """
    The results of an association test

    Attributes
    ----------
    data : npt.NDArray[np.float64]
        A numpy mixed array with fields: beta, pval, stderr

        It has shape num_variants x num_fields
    """

    data: npt.NDArray[np.float64, np.float64, np.float64]


class AssocTest(ABC):
    """
    Abstract class for performing phenotype-haplotype association tests

    Attributes
    ----------
    pval_thresh : float, optional
        The threshold of significance
    """
    def standardize(self, X: npt.NDArray[np.uint8]) -> npt.NDArray[np.float64]:
        """
        Standardize the genotypes so they have mean 0 and variance 1

        Parameters
        ----------
        X : npt.NDArray[np.float64]
            The genotypes, with shape n x p. There are only two dimensions.
            Each column is a haplotype and each row is a sample.

        Returns
        -------
        npt.NDArray[np.float64]
            An array with the same shape as X but standardized properly
        """
        std = np.std(X, axis=0)
        standardized = (X - np.mean(X, axis=0)) / std
        # # for variants where the stdev is 0, just set all values to 0 instead of nan
        # zero_elements = std == 0
        # standardized[:, zero_elements] = np.zeros(
        #     (X.shape[0], np.sum(zero_elements))
        # )
        return standardized

    @abstractmethod
    def run(
        self, X: npt.NDArray[np.float64], y: npt.NDArray[np.float64]
    ) -> AssocResults:
        """
        Run a series of phenotype-haplotype association tests for each haplotype
        (column) in X and return their p-values

        Parameters
        ----------
        X : npt.NDArray[np.float64]
            The genotypes, with shape n x p. There are only two dimensions.
            Each column is a haplotype and each row is a sample.
        y : npt.NDArray[np.float64]
            The phenotypes, with shape n x 1

        Returns
        -------
        npt.NDArray[np.float64]
            The p-values from testing each haplotype, with shape p x 1
        """
        pass


class AssocTestSimple(AssocTest):
    def __init__(self, with_bic=False):
        self.with_bic = with_bic
        super().__init__()

    def bic(self, n, residuals) -> float:
        """
        Return the BIC (Bayesian Information Criterion) for an OLS test

        This function follows the implementation in https://stackoverflow.com/a/58984868

        Parameters
        ----------
        n : int
            The number of observations
        residuals : npt.NDArray[np.float64]
            The residuals

        Returns
        -------
        float
            The Bayesian Information Criterion
        """
        k = 2
        sse = residuals.dot(residuals)
        if sse == 0:
            # if this is 0, then we already know that the BIC should be inf
            return float("inf")
        likelihood = -(n / 2) * (1 + np.log(2 * np.pi)) - (n / 2) * np.log(sse / n)
        return (-2 * likelihood) + ((k + 1) * np.log(n))

    def perform_test(
        self, x: npt.NDArray[np.float64], y: npt.NDArray[np.float64],
    ) -> tuple:
        """
        Perform the test for a single haplotype.

        Parameters
        ----------
        x : npt.NDArray[np.float64]
            The genotypes with shape n x 1 (for a single haplotype)
        y : npt.NDArray[np.float64]
            The phenotypes, with shape n x 1

        Returns
        -------
        tuple
            The slope, p-value, and stderr obtained from the test. The delta BIC is
            appended to the end if ``self.with_bic`` is True.
        """
        res = stats.linregress(x, y)
        if self.with_bic:
            y_hat = res.intercept + res.slope * x
            residuals = y - y_hat
            bic = self.bic(y.shape[0], residuals)
            return res.slope, res.pvalue, res.stderr, bic
        else:
            return res.slope, res.pvalue, res.stderr

    def run(
        self, X: npt.NDArray[np.float64], y: npt.NDArray[np.float64]
    ) -> AssocResults:
        """
        Implement AssocTest for a simple linear regression.

        Parameters
        ----------
        X : npt.NDArray[np.float64]
            The genotypes, with shape n x p. There are only two dimensions.
            Each row is a sample and each column is a haplotype.
        y : npt.NDArray[np.float64]
            The phenotypes, with shape n x 1

        Returns
        -------
        npt.NDArray[np.float64]
            The results from testing each haplotype, with shape p x 3
        """
        # use ordinary least squares for a simple regression
        # return an array of p-values
        return AssocResults(
            np.array(
                [
                    self.perform_test(X[:, variant_idx], y)
                    for variant_idx in range(X.shape[1])
                ],
                dtype=[
                    ("beta", np.float64),
                    ("pval", np.float64),
                    ("stderr", np.float64),
                ]
                + ([("bic", np.float64)] if self.with_bic else []),
            )
        )


class AssocTestSimpleSM(AssocTestSimple):
    def perform_test(
        self, x: npt.NDArray[np.float64], y: npt.NDArray[np.float64],
    ) -> tuple:
        """
        Perform the test for a single haplotype.

        Parameters
        ----------
        x : npt.NDArray[np.float64]
            The genotypes with shape n x 1 (for a single haplotype)
        y : npt.NDArray[np.float64]
            The phenotypes, with shape n x 1

        Returns
        -------
        tuple
            The slope, p-value, and stderr obtained from the test. The delta BIC is
            appended to the end if ``self.with_bic`` is True.
        """
        res = sm.OLS(y, sm.add_constant(x)).fit()
        params = res.params
        stderr = res.bse
        pvals = res.pvalues
        bic = res.bic
        if self.with_bic:
            return params[-1], pvals[-1], stderr[-1], bic
        else:
            return params[-1], pvals[-1], stderr[-1]


class AssocTestSimpleCovariates(AssocTestSimple):
    def __init__(self, covars: npt.NDArray[np.float64], with_bic=False):
        """
        Implement a subclass of AssocTestSimple with covariates.

        Parameters
        ----------
        covars : npt.NDArray[np.float64]
            Covariates to be included in the linear model. This should have shape
            n x m where each row (n) is sample and each column (m) is a covariate.
        """
        self.covars = covars
        super().__init__(with_bic=with_bic)

    def perform_test(self, x, y):
        """
        Perform the test for a single haplotype.

        Parameters
        ----------
        x : npt.NDArray[np.float64]
            The genotypes with shape n x 1 (for a single haplotype)
        y : npt.NDArray[np.float64]
            The phenotypes, with shape n x 1

        Returns
        -------
        tuple
            The slope, p-value, and stderr obtained from the test.
        """
        # the independent variables consist of this haplotype and the covariates
        X = np.hstack(x, self.covars)
        # initialize and create a multi-linear regression model
        mlr = sm.OLS(y, X)
        fit = mlr.fit()
        params = fit.params
        stderr = fit.bse
        pvals = fit.pvalues
        bic = fit.bic
        if self.with_bic:
            return params[0], pvals[0], stderr[0], bic
        else:
            return params[0], pvals[0], stderr[0]
