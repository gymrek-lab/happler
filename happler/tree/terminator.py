from __future__ import annotations
from logging import Logger
from abc import ABC, abstractmethod
from decimal import Decimal, getcontext

import numpy as np
from scipy.stats import t as t_dist
from haptools.logging import getLogger
from .corrector import Corrector, Bonferroni

from .assoc_test import NodeResults, NodeResultsExtra, AssocResults, AssocTest


class Terminator(ABC):
    """
    Abstract class for performing termination tests to check whether a branch should be
    terminated

    Attributes
    ----------
    thresh : float, optional
        The threshold of significance
    """

    def __init__(
        self, thresh: float = 0.05, corrector: Corrector = Bonferroni, log: Logger = None
    ):
        self.thresh = thresh
        self.corrector = corrector
        self.log = log or getLogger(self.__class__.__name__)
        super().__init__()

    @abstractmethod
    def compute_val(
        self,
        parent_res: NodeResults,
        node_res: NodeResults,
        results: AssocResults,
        best_idx: int,
        num_samps: int,
        num_tests: int,
        short_circuit: bool = True,
    ) -> tuple[float, float]:
        """
        Compute p-value or other thresholded value.
        This is a helper for the check() method

        Parameters
        ----------
        parent_res : NodeResults
            The results of the tests performed on the parent node
        node_res : NodeResults
            The results of the tests performed on the current node
        results: AssocResults
            All of the results for all of the variants tested on the current node
        best_idx : int
            The index of the best variant in the full set of results
        num_samps : int
            The number of samples tested
        num_tests : int
            The number of haplotypes tested
        short_circuit : bool, optional
            Should we produce values for all variants (F) or just the best ones (T)?

        Returns
        -------
        float
            The p-value
        float
            The t-score or BIC value
        """
        pass

    @abstractmethod
    def check(
        self,
        parent_res: NodeResults,
        node_res: NodeResults,
        results: AssocResults,
        best_idx: int,
        num_samps: int,
        num_tests: int,
        logger: Logger,
    ) -> bool:
        """
        Check whether a branch should be terminated.

        Parameters
        ----------
        parent_res : NodeResults
            The results of the tests performed on the parent node
        node_res : NodeResults
            The results of the tests performed on the current node
        results: AssocResults
            All of the results for all of the variants tested on the current node
        best_idx : int
            The index of the best variant in the full set of results
        num_samps : int
            The number of samples tested
        num_tests : int
            The number of haplotypes tested
        log: Logger
            A logging instance for recording debug statements.

        Returns
        -------
        bool
            True if the branch should be terminated, False otherwise
        """
        pass


class TTestTerminator(Terminator):
    def compute_val(
        self,
        parent_res: NodeResults,
        node_res: NodeResults,
        results: AssocResults,
        best_idx: int,
        num_samps: int,
        num_tests: int,
        parent_corr: float = 0,
        short_circuit: bool = True,
    ) -> tuple[float, float]:
        t_stat = None
        if parent_res:
            # before we do any calculations, check whether the effect sizes have
            # improved and return True if they haven't
            if short_circuit and (
                np.isnan(node_res.beta)
                or (np.abs(node_res.beta) - np.abs(parent_res.beta)) <= 0
            ):
                # terminate if the effect sizes have gone in the opposite direction
                self.log.debug("Terminated b/c effect size did not improve")
                return True
            # perform a one tailed, two-sample t-test using the difference of the effect sizes
            # first, we compute the covariance of the effect sizes
            cov = 0
            if parent_corr is not None:
                cov = parent_corr * parent_res.stderr * results.data["stderr"]
            # now, we can compute the standard error of the difference of the effect sizes
            std_err = ((results.data["stderr"] ** 2) + (parent_res.stderr**2)) - 2 * cov
            if std_err[best_idx] == 0:
                # if we have a standard error less than 0, then we already know the result is
                # significant! It doesn't matter what the effect sizes are b/c t_stat
                # will be inf
                return (0, np.inf)
            std_err = np.sqrt(std_err / 2)
            # then, we compute the test statistic
            # use np.abs to account for the directions that the effect size may take
            t_stat = (np.abs(results.data["beta"]) - np.abs(parent_res.beta)) / std_err
            # use a one-tailed test here b/c either the effect size becomes more
            # negative or it becomes more positive
            dof = 2 * (num_samps - 2)
            pval = t_dist.sf(t_stat, df=dof)
            if (pval == 0).any():
                pval = pval.astype(object)
                # retrieve the pvals at a higher precision as Decimal objects
                pval[pval == 0] = [
                    AssocTest.pval_as_decimal(tval, dof, precision=10)
                    for tval in t_stat[pval == 0]
                ]
        else:
            # parent_res = None when the parent node is the root node
            pval = results.data["pval"]
        if self.corrector is not None:
            corrector = self.corrector(thresh=self.thresh, log=self.log)
            pval = corrector.correct(pval, num_samps, len(results.data))[best_idx]
        else:
            pval = pval[best_idx]
        if t_stat is not None:
            t_stat = t_stat[best_idx]
        return (pval, t_stat)

    def check(
        self,
        parent_res: NodeResults,
        node_res: NodeResults,
        results: AssocResults,
        best_idx: int,
        num_samps: int,
        num_tests: int,
        parent_corr: float = 0,
    ) -> bool:
        computed_val = self.compute_val(
            parent_res,
            node_res,
            results,
            best_idx,
            num_samps,
            num_tests,
            parent_corr,
        )
        if isinstance(computed_val, bool):
            return computed_val
        else:
            pval, t_stat = computed_val
        if not isinstance(pval, Decimal) and np.isnan(pval):
            raise ValueError(
                "Encountered an nan p-value! Check your data for irregularities."
            )
        if pval >= self.thresh:
            self.log.debug(
                f"Terminated with t-stat {t_stat} and p-value {pval} >= {self.thresh}"
            )
            return True
        self.log.debug(
            f"Significant with t-stat {t_stat} and p-value {pval} < {self.thresh}"
        )
        return False


class BICTerminator(Terminator):
    def __init__(self, thresh: float = 0.05, bic_thresh: float = 3, log: Logger = None):
        super().__init__()
        self.thresh = thresh
        self.bic_thresh = bic_thresh
        self.log = log or getLogger(self.__class__.__name__)

    def compute_val(
        self,
        parent_res: NodeResultsExtra,
        node_res: NodeResultsExtra,
        results: AssocResults,
        best_idx: int,
        num_samps: int,
        num_tests: int,
        parent_corr: float = 0,
        short_circuit: bool = True,
    ) -> tuple[float, float]:
        stat = None
        pval = None
        if parent_res:
            # before we do any calculations, check whether the effect sizes have
            # improved and return True if they haven't
            if short_circuit and (
                np.isnan(node_res.beta)
                or (np.abs(node_res.beta) - np.abs(parent_res.beta)) <= 0
            ):
                # terminate if the effect sizes have gone in the opposite direction
                self.log.debug("Terminated b/c effect size did not improve")
                return True
            stat = parent_res.bic - results.data["bic"]
            # compute the bayes factor approximation from the delta BIC:
            # https://easystats.github.io/bayestestR/reference/bic_to_bf.html
            stat = np.exp(stat / 2)
            stat = stat[best_idx]
        else:
            # parent_res = None when the parent node is the root node
            pval = results.data["pval"]
            # correct for multiple hypothesis testing
            if self.corrector is not None:
                # I dunno why, but this needs None as the first arg for some reason
                pval = self.corrector.correct(None, pval, num_samps, len(pval))[best_idx]
            else:
                pval = pval[best_idx]
        return pval, stat

    def check(
        self,
        parent_res: NodeResultsExtra,
        node_res: NodeResultsExtra,
        results: AssocResults,
        best_idx: int,
        num_samps: int,
        num_tests: int,
    ):
        computed_val = self.compute_val(
            parent_res,
            node_res,
            results,
            best_idx,
            num_samps,
            num_tests,
        )
        if isinstance(computed_val, bool):
            return computed_val
        else:
            pval, stat = computed_val
        if stat is None:
            # TODO: handle this case by using delta BIC to rank, instead?
            if pval >= self.thresh:
                self.log.debug(
                    f"Terminated with delta BIC {stat} and p-value {pval} >= {self.thresh}"
                )
                return True
        else:
            if stat < self.bic_thresh:
                self.log.debug(
                    f"Terminated with delta BIC {stat} < {self.thresh} and p-value {pval}"
                )
                return True
        self.log.debug(f"Significant with delta BIC {stat} > {self.thresh}")
        return False
