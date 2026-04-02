from __future__ import annotations
from logging import getLogger
from dataclasses import dataclass
from abc import ABC, abstractmethod
from decimal import Decimal, getcontext

import numpy as np
from scipy import stats
import numpy.typing as npt
import statsmodels.api as sm

# Import JAX for JIT compilation
import jax
import jax.numpy as jnp

# Configure JAX for 64-bit precision to match NumPy behavior
jax.config.update("jax_enable_x64", True)


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


@dataclass(frozen=True)
class NodeResultsBIC:
    """
    The results of testing SNPs at a node in the tree

    Attributes
    ----------
    bic : float
        The best BIC among all of the SNPs
    """

    bic: float

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

    data: npt.NDArray


class AssocTest(ABC):
    """
    Abstract class for performing phenotype-haplotype association tests

    Attributes
    ----------
    with_bic: bool
        Whether to also output the BIC (Bayesian Information Criteria)
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
        # for variants where the stdev is 0, just set all values to 0 instead of nan
        zero_elements = std == 0
        standardized[:, zero_elements] = np.zeros((X.shape[0], np.sum(zero_elements)))
        return standardized

    @staticmethod
    def pval_as_decimal(t_stat: float, df: int, precision: int = 1000) -> Decimal:
        """
        Given a t statistic, return the associated p-value as a high precision Decimal

        This can be helpful when a p-value is too small to be represented by the
        precision in python float64s

        Parameters
        ----------
        t_stat : float
            The t statistic from a simple linear model
        df : int
            The degrees of freedom of the associated t distribution
        precision: int
            The precision of the returned value

        Raises
        ------
        ValueError
            This function is only valid when the t distribution is approximately
            equivalent to the normal distribution. Thus, if df < 1000, we raise a
            ValueError to indicate that the function will not return an accurate value

        Returns
        -------
        Decimal
            An approximate, higher precision p-value for the provided t statistic
        """
        if df < 1000:
            getLogger().warning(
                "You need a larger sample size to approximate this p-value"
            )
        log10_pval = stats.norm.logsf(np.abs(t_stat)) / np.log(10) + np.log10(2)
        # set the desired precision
        getcontext().prec = precision
        # calculate p-value
        return Decimal(10) ** Decimal(log10_pval)

    @abstractmethod
    def run(self, X: npt.NDArray[np.float64], y: npt.NDArray[np.float64]) -> AssocResults:
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
            ("pval", object),
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
        raise ValueError(
            "BIC is not properly implemented in this class. "
            "It will yield weird/incorrect values."
        )
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
            The slope, p-value, and stderr obtained from the test. The BIC is
            appended to the end if ``self.with_bic`` is True.
        """
        try:
            res = stats.linregress(x, y)
        except ValueError:
            # This can happen, for example, when all x values are identical
            # In this case, we just output the worst possible values
            if self.with_bic:
                return 0, 1, np.inf, 0
            return 0, 1, np.inf
        if self.with_bic:
            y_hat = res.intercept + res.slope * x
            residuals = y - y_hat
            bic = self.bic(y.shape[0], residuals)
            return res.slope, res.pvalue, res.stderr, bic
        else:
            return res.slope, res.pvalue, res.stderr

    def run(self, X: npt.NDArray[np.float64], y: npt.NDArray[np.float64]) -> AssocResults:
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
            The slope, p-value, and stderr obtained from the test. The BIC is
            appended to the end if ``self.with_bic`` is True.
        """
        res = self._get_sm_result(x, y)
        param = res.params[-1]
        stderr = res.bse[-1]
        pval = res.pvalues[-1]
        bic = res.bic
        # handle cases where our p-values are too powerful
        if pval == 0:
            # retrieve the pval at a higher precision
            pval = AssocTestSimpleSM.pval_as_decimal(
                res.tvalues[-1], res.df_resid, precision=10
            )
        if self.with_bic:
            return param, pval, stderr, bic
        else:
            return param, pval, stderr


@jax.jit
def _compute_bic_jit(X: jnp.ndarray, yc: jnp.ndarray) -> jnp.ndarray:
    """
    JIT-compiled helper function for computing BIC in vectorized OLS regression.

    This function uses JAX for JIT compilation and GPU/TPU acceleration when available.
    Includes explicit inf handling to ensure numerical equivalence with NumPy, as JAX
    JIT optimization can sometimes change exact numerics (e.g., log(exp(x)) -> x).

    Parameters
    ----------
    X : jnp.ndarray
        The genotypes with shape (n, p) where n is samples and p is variants
    yc : jnp.ndarray
        The centered phenotypes with shape (n, 1)

    Returns
    -------
    jnp.ndarray
        The BIC values with shape (p,)
    """
    n = X.shape[0]
    nobs2 = n / 2.0
    log2pi = jnp.log(2 * jnp.pi)

    # Center X
    xc = X - jnp.mean(X, axis=0)  # (n, p)

    # Vectorized simple OLS with intercept
    sxx = jnp.sum(xc**2, axis=0)  # (p,)
    sxy = jnp.sum(xc * yc, axis=0)  # (p,)

    # Compute slopes, handling division by zero
    b1 = jnp.where(sxx > 0, sxy / sxx, 0.0)  # slopes, (p,)

    # Residuals for each column's model
    syy = jnp.sum(yc**2)  # scalar
    ssr = syy - 2 * b1 * sxy + (b1**2) * sxx

    # Explicitly handle cases that should produce inf values
    # When SSR is very small or zero, log(SSR/n) should be -inf
    # JAX JIT might optimize this differently than NumPy
    ssr_threshold = 1e-300  # Below this, treat as zero
    ssr_safe = jnp.where(ssr > ssr_threshold, ssr, ssr_threshold)

    # statsmodels-style profile log-likelihood per column
    # ll = -n/2 * [ log(2π) + log(SSR/n) + 1 ]
    ll = -nobs2 * (log2pi + jnp.log(ssr_safe / n) + 1.0)  # (p,)

    # Number of parameters k: intercept + slope = 2
    bic = -2 * ll + 2 * jnp.log(n)  # (p,)

    # Explicitly set to -inf where SSR was effectively zero
    bic = jnp.where(ssr <= ssr_threshold, jnp.array(-jnp.inf), bic)

    return bic


class AssocTestSimpleFastBIC(AssocTestSimpleSM):
    """
    Calculate only BIC in a quick, vectorized fashion without statsmodels
    """

    def __init__(self, chunk_size: int = None):
        """
        Override the parent's __init__
        """
        self.results_type = NodeResultsBIC
        self.chunk_size = chunk_size

    def _compute_bic(X: npt.NDArray[np.float64], yc: npt.NDArray[np.float64]) -> npt.NDArray:
        n = X.shape[0]
        nobs2 = n / 2.0
        log2pi = np.log(2 * np.pi)

        # Center X and y
        xc = X - X.mean(axis=0)  # (n, p)

        # Vectorized simple OLS with intercept
        sxx = np.sum(xc**2, axis=0)  # (p,)
        sxy = np.sum(xc * yc, axis=0)  # (p,)

        with np.errstate(divide="ignore", invalid="ignore"):
            b1 = np.where(sxx > 0, sxy / sxx, 0.0)  # slopes, (p,)

        # Residuals for each column's model: r = (y - ym) - b1 * (X - xm)
        syy = float(np.sum(yc**2))  # scalar
        ssr = syy - 2 * b1 * sxy + (b1**2) * sxx

        # statsmodels-style profile log-likelihood per column
        # ll = -n/2 * [ log(2π) + log(SSR/n) + 1 ]
        with np.errstate(divide="ignore"):
            ll = -nobs2 * (log2pi + np.log(ssr / n) + 1.0)  # (p,)

        # Number of parameters k: intercept + slope = 2
        return -2 * ll + 2 * np.log(n)  # (p,)

    def perform_test(
        self, X: npt.NDArray[np.float64], yc: npt.NDArray[np.float64]
    ) -> npt.NDArray:
        """
        Perform the test for a chunk of haplotypes using JAX JIT compilation

        This method uses the JIT-compiled _compute_bic_jit helper function
        for improved computational efficiency with GPU/TPU acceleration.

        Parameters
        ----------
        X : npt.NDArray[np.float64]
            The genotypes with shape (n, p)
        yc : npt.NDArray[np.float64]
            The phenotypes, with shape (n, 1)
            They are assumed to be centered already

        Returns
        -------
        npt.NDArray[np.float64]
            The BIC values from testing this chunk of haplotypes, with shape (p,)
        """
        use_jax = False
        if use_jax:
            # Call JAX JIT-compiled helper and convert result to writable NumPy array
            return np.array(_compute_bic_jit(X, yc), copy=True)
        else:
            return self._compute_bic(X, yc)

    def run(self, X: npt.NDArray[np.float64], y: npt.NDArray[np.float64]) -> AssocResults:
        """
        Implement AssocTest for a simple, univariate OLS model: y ~ 1 + X[:, j]

        Does not use statsmodels at all but replicates its behavior

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
            The results from testing each haplotype, with shape p x 1
        """
        if len(y.shape) != 2:
            y = y[:, np.newaxis]  # (n, 1)

        yc = y - float(y.mean())  # (n, 1)
        bic_vals = np.zeros((X.shape[1]), dtype=np.float64)  # (p, 1)

        chunks = self.chunk_size
        if chunks is None or chunks > len(bic_vals):
            chunks = len(bic_vals)

        for start in range(0, len(bic_vals), chunks):
            end = start + chunks
            if end > len(bic_vals):
                end = len(bic_vals)
            size = end - start

            bic_vals[start:end] = self.perform_test(X[:, start:end], yc)

        return AssocResults(bic_vals.astype([("bic", np.float64)]))


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
        parent_corr: float = 0,
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
            The slope, p-value, and stderr obtained from the test. The BIC is
            appended to the end if ``self.with_bic`` is True.
        """
        if self.with_bic:
            beta, pval, stderr, bic = super().perform_test(x, y)
        else:
            beta, pval, stderr = super().perform_test(x, y)
        if parent_res is None:
            t_score = 0
        else:
            cov = parent_corr * parent_res.stderr * stderr
            if np.isnan(cov):
                cov = 0
            std_err = np.sqrt((((stderr**2) + (parent_res.stderr**2)) / 2) - 2 * cov)
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
        parent_corr: npt.NDArray[np.float64] = None,
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
        if parent_corr is None:
            corr = lambda i: 0
        else:
            corr = lambda i: parent_corr[i]
        # use ordinary least squares for a simple regression
        # return an array of p-values
        return AssocResults(
            np.array(
                [
                    self.perform_test(X[:, variant_idx], y, parent_res, corr(variant_idx))
                    for variant_idx in range(X.shape[1])
                ],
                dtype=self.return_dtype,
            )
        )
