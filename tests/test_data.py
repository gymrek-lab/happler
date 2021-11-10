import numpy as np
from pathlib import Path
from pytest import approx

from happler.data import Genotypes


DATADIR = Path(__file__).parent.joinpath('data')

def test_load_genotypes():
    expected = np.zeros(60).reshape((4,5,3))
    expected[:, :, 2] = 1

    gts = Genotypes(DATADIR.joinpath('simple.vcf'))
    gts.load()
    np.testing.assert_allclose(gts.data, expected)

    gts.check_phase()
    expected = expected[:, :, :2]
    np.testing.assert_allclose(gts.data, expected)

    gts.to_MAC()
    np.testing.assert_allclose(gts.data, expected)
