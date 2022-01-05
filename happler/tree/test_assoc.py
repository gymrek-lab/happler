from __future__ import annotations
from pathlib import Path
from abc import ABC, abstractmethod

import numpy as np
from scipy import stats
import numpy.typing as npt


class TestAssoc(ABC):
    """
    Abstract class for performing phenotype-haplotype association tests

    Attributes
    ----------
    pval_thresh : float, optional
        The threshold of significance
    """

    def __init__(self, pval_thresh: float = 0.05):
        self.pval_thresh = pval_thresh
        super().__init__()

    @abstractmethod
    def run(
        self, X: npt.NDArray[np.float64], y: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
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


class TestAssocSimple(TestAssoc):
    def run(
        self, X: npt.NDArray[np.float64], y: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """
        Implement TestAssoc for a simple linear regression.

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
        # use ordinary least squares for a simple regression
        # return an array of p-values
        return np.array(
            [
                stats.linregress(X[:, variant_idx], y).pvalue
                for variant_idx in range(X.shape[1])
            ],
            dtype=np.float64,
        )
