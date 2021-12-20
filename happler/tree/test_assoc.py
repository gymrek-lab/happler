from __future__ import annotations
from pathlib import Path
from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt
import statsmodels.api as sm


class TestAssoc(ABC):

    @abstractmethod
    def run(self, X: npt.NDArray[np.float64], y: npt.NDArray[np.float64]):
        pass


class TestAssocSimple(TestAssoc):

    def run(self, X: npt.NDArray[np.float64], y: npt.NDArray[np.float64]):
        # TODO: the current code fits a multivariate regression model
        # but we want to fit each haplotype separately
        mod = sm.OLS(y, X)
        res = mod.fit()
        p_values = res.summary2().tables[1]['P>|t|']
        return p_values
