from __future__ import annotations
from logging import Logger
from abc import ABC, abstractmethod

import numpy as np
from scipy.stats import t as t_dist
from haptools.logging import getLogger
from .corrector import Corrector, Bonferroni

from .assoc_test import NodeResults, NodeResultsExtra, AssocResults


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
    def check(
        self,
        parent_res: NodeResults,
        node_res: NodeResults,
        results: AssocResults,
        best_idx: int,
        num_samps: int,
        num_tests: int,
    ) -> bool:
        t_stat = None
        if parent_res:
            # before we do any calculations, check whether the effect sizes have
            # improved and return True if they haven't
            if (
                np.isnan(node_res.beta)
                or (np.abs(node_res.beta) - np.abs(parent_res.beta)) <= 0
            ):
                # terminate if the effect sizes have gone in the opposite direction
                self.log.debug("Terminated b/c effect size did not improve")
                return True
            # perform a two tailed, two-sample t-test using the difference of the effect sizes
            # first, we compute the standard error of the difference of the effect sizes
            std_err = np.sqrt(
                ((results.data["stderr"] ** 2) + (parent_res.stderr**2)) / 2
            )
            if std_err[best_idx] == 0:
                # if we have a standard error of 0, then we already know the result is
                # significant! It doesn't matter what the effect sizes are b/c t_stat
                # will be inf
                return False
            # then, we compute the test statistic
            # use np.abs to account for the directions that the effect size may take
            t_stat = (np.abs(results.data["beta"]) - np.abs(parent_res.beta)) / std_err
            # use a one-tailed test here b/c either the effect size becomes more
            # negative or it becomes more positive
            pval = t_dist.sf(t_stat, df=2 * (num_samps - 2))
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
        if np.isnan(pval):
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
    def __init__(self, thresh: float = 10, log: Logger = None):
        super().__init__()
        self.thresh = thresh
        self.log = log or getLogger(self.__class__.__name__)

    def check(
        self,
        parent_res: NodeResultsExtra,
        node_res: NodeResultsExtra,
        results: AssocResults,
        best_idx: int,
        num_samps: int,
        num_tests: int,
    ) -> bool:
        stat = None
        if parent_res:
            # before we do any calculations, check whether the effect sizes have
            # improved and return True if they haven't
            if (
                np.isnan(node_res.beta)
                or (np.abs(node_res.beta) - np.abs(parent_res.beta)) <= 0
            ):
                # terminate if the effect sizes have gone in the opposite direction
                self.log.debug("Terminated b/c effect size did not improve")
                return True
            stat = results.data["bic"] - parent_res.bic
            # just choose an arbitrary threshold
            if stat[best_idx] < self.thresh:
                self.log.debug("Terminated with delta BIC {}".format(stat))
                return True
        else:
            # parent_res = None when the parent node is the root node
            pval = results.data["pval"]
            # correct for multiple hypothesis testing
            if self.corrector is not None:
                # I dunno why, but this needs None as the first arg for some reason
                pval = self.corrector.correct(None, pval, num_samps, len(pval))[best_idx]
            else:
                pval = pval[best_idx]
            if pval >= self.thresh:
                self.log.debug("Terminated with p-value {}".format(pval))
                return True
        self.log.debug(
            "Significant with "
            + ("delta BIC {}".format(stat) if parent_res else "pval {}".format(pval))
        )
        return False
