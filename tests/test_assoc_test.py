from pathlib import Path

import pytest
import numpy as np
from haptools.data import GenotypesPLINK, Phenotypes

from happler.tree import (
    AssocTestSimpleSM,
    AssocTestSimpleFastBIC,
    NodeResultsBIC,
    NodeResultsExtra,
)

DATADIR = Path(__file__).parent.joinpath("data")


def test_bic_methods():
    """
    Test that our custom method of computing BIC is the same as statsmodels's
    """
    for dataset in (
        "15_38177344-39177344_1000G",
        "19_45401409-46401409_1000G",
        "19_45401409-46401409_1000G.multi",
    ):
        gt = GenotypesPLINK.load(DATADIR / (dataset + ".pgen")).data.sum(axis=2)
        pheno = Phenotypes.load(DATADIR / (dataset + ".pheno"))
        pt = pheno.data

        tester = AssocTestSimpleSM(with_bic=True)
        res_sm = tester.run(gt, pt)
        tester = AssocTestSimpleFastBIC()
        res_fast = tester.run(gt, pt)

        np.testing.assert_allclose(res_sm.data["bic"], res_fast.data["bic"])

        # also test after standardizing
        pheno.standardize()
        pt = pheno.data

        res_fast = tester.run(gt, pt)
        np.testing.assert_allclose(res_sm.data["bic"], res_fast.data["bic"])

        # and what if we also standardize the genotypes?
        gt = tester.standardize(gt)

        res_fast = tester.run(gt, pt)
        np.testing.assert_allclose(res_sm.data["bic"], res_fast.data["bic"])
