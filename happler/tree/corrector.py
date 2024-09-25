from __future__ import annotations
from logging import Logger
from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt
from haptools.logging import getLogger
from statsmodels.stats.multitest import fdrcorrection

# try to import scipy.stats.false_discovery_control
# which is only available when py >= 3.9 and scipy >= 1.11.0
try:
    from scipy.stats import false_discovery_control
except ImportError:
    false_discovery_control = None


class Corrector(ABC):
    """
    Abstract class for correcting p-values due to multiple tests

    Attributes
    ----------
    thresh : float, optional
        The threshold of significance
    """

    def __init__(self, thresh: float = 0.05, log: Logger = None):
        self.thresh = thresh
        self.log = log or getLogger(self.__class__.__name__)
        super().__init__()

    @abstractmethod
    def correct(
        self,
        pvals: npt.NDArray,
        num_samps: int,
        num_tests: int,
    ) -> npt.NDArray:
        """
        Correct a set of p-values

        Parameters
        ----------
        pvals: npt.NDArray
            The p-values to be corrected
        num_samps : int
            The number of samples tested
        num_tests : int
            The number of haplotypes tested
        log: Logger
            A logging instance for recording debug statements.

        Returns
        -------
        npt.NDArray
            A set of corrected p-values in the same order as the input array
        """
        pass


class Bonferroni(Corrector):
    def correct(
        self,
        pvals: npt.NDArray,
        num_samps: int,
        num_tests: int,
    ) -> npt.NDArray:
        """
        Refer to the documentation of :py:meth:`~.Corrector.correct`

        This class implements the Bonferroni correction
        """
        pvals = pvals * num_tests
        pvals[pvals > 1] = 1
        return pvals


class BH(Corrector):
    if false_discovery_control is None:

        def correct(
            self,
            pvals: npt.NDArray,
            num_samps: int,
            num_tests: int,
        ) -> npt.NDArray:
            """
            Refer to the documentation of :py:meth:`~.Corrector.correct`
            """
            raise NotImplementedError

    else:

        def correct(
            self,
            pvals: npt.NDArray,
            num_samps: int,
            num_tests: int,
        ) -> npt.NDArray:
            """
            Refer to the documentation of :py:meth:`~.Corrector.correct`

            This class implements the Benjamini-Hochberg procedure
            """
            return false_discovery_control(pvals, method="bh")


class BHSM(BH):
    def correct(
        self,
        pvals: npt.NDArray,
        num_samps: int,
        num_tests: int,
    ) -> npt.NDArray:
        """
        Refer to the documentation of :py:meth:`~.Corrector.correct`

        This class implements the Benjamini-Hochberg procedure with statsmodels
        """
        return fdrcorrection(pvals, alpha=self.thresh, method="poscorr")[1]
