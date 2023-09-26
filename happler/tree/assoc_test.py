from __future__ import annotations
from dataclasses import dataclass
from abc import ABC, abstractmethod

import numpy as np
from scipy import stats
import numpy.typing as npt
import statsmodels.api as sm
from scipy.stats import t as t_dist


# We declare this class to be a dataclass to automatically define __init__ and a few
# other methods. We use frozen=True to make it immutable.
@dataclass(frozen=True)
class NodeResults:
    """
    The results of testing SNPs at a node in the tree

    Attributes
    ----------
    beta : float
        The best effect size among all of the SNPs tried
    pval : float
        The best p-value among all of the SNPs tried
    stderr: float
        The standard error of beta
    """

    beta: float
    pval: float
    stderr: float

    def __getitem__(self, item):
        """
        Define a getter so that we can access elements like this:

        ``obj['field_name']``

        in addition to this:

        ``obj.field_name``
        """
        return getattr(self, item)

    def __repr__(self):
        return (
            "{" + ", ".join("{}={:.2e}".format(*i) for i in self.__dict__.items()) + "}"
        )

    @classmethod
    def from_np(cls, np_mixed_arr_var: np.void) -> NodeResults:
        class_attributes = cls.__dict__["__dataclass_fields__"].keys()
        return cls(**dict(zip(class_attributes, np_mixed_arr_var)))


@dataclass(frozen=True, repr=False)
class NodeResultsExtra(NodeResults):
    bic: float


@dataclass(frozen=True, repr=False)
class NodeResultsTScore(NodeResults):
    tscore: float


@dataclass(frozen=True, repr=False)
class NodeResultsExtraTScore(NodeResultsExtra):
    tscore: float


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
    def __init__(self, with_bic: bool = False):
        self.with_bic = with_bic
        self.results_type = NodeResults
        if self.with_bic:
            self.results_type = NodeResultsExtra

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
        zero_elements = std == 0
        standardized[:, zero_elements] = np.zeros((X.shape[0], np.sum(zero_elements)))
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
    def __init__(self, with_bic: bool = False):
        super().__init__(with_bic=with_bic)
        self.return_dtype = [
            ("beta", np.float64),
            ("pval", np.float64),
            ("stderr", np.float64),
        ]
        if self.with_bic:
            self.return_dtype.append(("bic", np.float64))

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
        self,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
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
        X = self.standardize(X)
        # use ordinary least squares for a simple regression
        # return an array of p-values
        return AssocResults(
            np.array(
                [
                    self.perform_test(X[:, variant_idx], y)
                    for variant_idx in range(X.shape[1])
                ],
                dtype=self.return_dtype,
            )
        )


class AssocTestSimpleSM(AssocTestSimple):
    def _get_sm_result(
        self,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
    ) -> sm.regression.linear_model.RegressionResults:
        return sm.OLS(y, sm.add_constant(x)).fit()

    def perform_test(
        self,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
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
        res = self._get_sm_result(x, y)
        params = res.params
        stderr = res.bse
        pvals = res.pvalues
        bic = res.bic
        if self.with_bic:
            return params[-1], pvals[-1], stderr[-1], bic
        else:
            return params[-1], pvals[-1], stderr[-1]


class AssocTestSimpleCovariates(AssocTestSimpleSM):
    def __init__(self, covars: npt.NDArray[np.float64], with_bic: bool = False):
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

    def _get_sm_result(
        self,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
    ) -> sm.regression.linear_model.RegressionResults:
        X = sm.add_constant(np.column_stack((self.covars, x)))
        return sm.OLS(y, X).fit()


class AssocTestSimpleSMTScore(AssocTestSimpleSM):
    def __init__(self, with_bic: bool = False):
        """
        Implement a subclass of AssocTestSimple with covariates.

        Parameters
        ----------
        covars : npt.NDArray[np.float64]
            Covariates to be included in the linear model. This should have shape
            n x m where each row (n) is sample and each column (m) is a covariate.
        """
        super().__init__(with_bic=with_bic)
        self.return_dtype.append(("tscore", np.float64))
        self.results_type = NodeResultsTScore
        if self.with_bic:
            self.results_type = NodeResultsExtraTScore

    def perform_test(
        self,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        parent_res: NodeResults = None,
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
        if self.with_bic:
            beta, pval, stderr, bic = super().perform_test(x, y)
        else:
            beta, pval, stderr = super().perform_test(x, y)
        if parent_res is None:
            t_score = 0
        else:
            std_err = np.sqrt(((stderr**2) + (parent_res.stderr**2)) / 2)
            t_score = (np.abs(beta) - np.abs(parent_res.beta)) / std_err
        if self.with_bic:
            return beta, pval, stderr, bic, t_score
        else:
            return beta, pval, stderr, t_score

    def run(
        self,
        X: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        parent_res: NodeResults = None,
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
        X = self.standardize(X)
        # use ordinary least squares for a simple regression
        # return an array of p-values
        return AssocResults(
            np.array(
                [
                    self.perform_test(X[:, variant_idx], y, parent_res)
                    for variant_idx in range(X.shape[1])
                ],
                dtype=self.return_dtype,
            )
        )
