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
        A numpy mixed array with fields: beta and pval
    """

    data: npt.NDArray[np.float64, np.float64]


class AssocTest(ABC):
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
            The p-values from testing each haplotype, with shape p x 1
        """
        extract_vals = lambda linreg: (linreg.slope, linreg.pvalue)
        # use ordinary least squares for a simple regression
        # return an array of p-values
        return AssocResults(
            np.array(
                [
                    extract_vals(stats.linregress(X[:, variant_idx], y))
                    for variant_idx in range(X.shape[1])
                ],
                dtype=[
                    ("beta", np.float64),
                    ("pval", np.float64),
                ],
            )
        )
