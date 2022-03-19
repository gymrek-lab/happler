from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass
from abc import ABC, abstractmethod

import numpy as np
from scipy import stats
import numpy.typing as npt


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

        This function follows the implementatin in https://stackoverflow.com/a/58984868

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
            return float('inf')
        likelihood = -(n / 2) * (1 + np.log(2 * np.pi)) - (n / 2) * np.log(sse / n)
        return (-2 * likelihood) + ((k + 1) * np.log(n))

    def perform_test(self, X, y):
        res = stats.linregress(X, y)
        if self.with_bic:
            y_hat = res.intercept + res.slope * X
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
            Each column is a haplotype and each row is a sample.
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
