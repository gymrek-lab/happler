from __future__ import annotations
from logging import Logger
from abc import ABC, abstractmethod

import numpy as np
from scipy.stats import t as t_dist

from .tree import NodeResults, NodeResultsExtra


class Terminator(ABC):
    """
    Abstract class for performing termination tests to check whether a branch should be
    terminated

    Attributes
    ----------
    thresh : float, optional
        The threshold of significance
    """

    def __init__(self, thresh: float = 0.05):
        self.thresh = thresh
        super().__init__()

    @abstractmethod
    def check(
        self,
        parent_res: NodeResults,
        node_res: NodeResults,
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
        num_samps: int,
        num_tests: int,
        logger: Logger,
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
                logger.debug("Terminated b/c effect size went in wrong direction")
                return True
            # perform a two tailed, two-sample t-test using the difference of the effect sizes
            # first, we compute the standard error of the difference of the effect sizes
            std_err = np.sqrt(((node_res.stderr**2) + (parent_res.stderr**2)) / 2)
            if std_err == 0:
                # if we have a standard error of 0, then we already know the result is
                # significant! It doesn't matter what the effect sizes are b/c t_stat
                # will be inf
                return False
            # then, we compute the test statistic
            # use np.abs to account for the directions that the effect size may take
            t_stat = np.abs(node_res.beta - parent_res.beta) / std_err
            # use a one-tailed test here b/c either the effect size becomes more
            # negative or it becomes more positive
            pval = t_dist.cdf(-t_stat, df=2 * (num_samps - 2))
        else:
            # parent_res = None when the parent node is the root node
            pval = node_res.pval
        if np.isnan(pval):
            raise ValueError(
                "Encountered an nan p-value! Check your data for irregularities."
            )
        # correct for multiple hypothesis testing
        # For now, we use the Bonferroni correction
        if pval >= (self.thresh / num_tests):
            logger.debug(
                "Terminated with t-stat {} and p-value {}".format(t_stat, pval)
            )
            return True
        logger.debug("Significant with t-stat {} and p-value {}".format(t_stat, pval))
        return False


class BICTerminator(Terminator):
    def __init__(self, thresh: float = 10):
        super().__init__()
        self.thresh = thresh

    def check(
        self,
        parent_res: NodeResultsExtra,
        node_res: NodeResultsExtra,
        num_samps: int,
        num_tests: int,
        logger: Logger,
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
                logger.debug("Terminated b/c effect size went in wrong direction")
                return True
            stat = node_res.bic - parent_res.bic
            # just choose an arbitrary threshold
            if stat < self.thresh:
                logger.debug("Terminated with delta BIC {}".format(stat))
                return True
        else:
            # parent_res = None when the parent node is the root node
            pval = node_res.pval
            # correct for multiple hypothesis testing
            # For now, we use the Bonferroni correction
            if pval >= (self.thresh / num_tests):
                logger.debug("Terminated with p-value {}".format(pval))
                return True
        logger.debug(
            "Significant with "
            + ("delta BIC {}".format(stat) if parent_res else "pval {}".format(pval))
        )
        return False
