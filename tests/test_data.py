import pytest
import numpy as np
from pathlib import Path

from happler.data import Genotypes, Phenotypes


DATADIR = Path(__file__).parent.joinpath("data")


def test_load_genotypes():
    # create a GT matrix with shape: samples x SNPs x (strands+phase)
    expected = np.zeros(60).reshape((5, 4, 3)).astype(np.bool_)
    expected[:4, 1, 1] = 1
    expected[2:4, 1, 0] = 1
    expected[:, :, 2] = 1

    # can we load the data from the VCF?
    gts = Genotypes(DATADIR.joinpath("simple.vcf"))
    gts.load()
    np.testing.assert_allclose(gts.data, expected)
    assert gts.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

    # try loading the data again - it should fail b/c we've already done it
    with pytest.raises(AssertionError):
        gts.load()

    # force one of the het SNPs to be unphased and check that we get an error message
    gts.data[1, 1, 2] = 0
    with pytest.raises(ValueError) as info:
        gts.check_phase()
    assert (
        str(info.value)
        == "Variant with ID 1:10116:A:G at POS 1:10116 is unphased for sample HG00097"
    )
    gts.data[1, 1, 2] = 1

    # check phase and remove the phase axis
    gts.check_phase()
    expected = expected[:, :, :2]
    np.testing.assert_allclose(gts.data, expected)

    # try to check phase again - it should fail b/c we've already done it before
    with pytest.raises(AssertionError):
        gts.check_phase()

    # convert the matrix of alt allele counts to a matrix of minor allele counts
    assert gts.variants["aaf"][1] == 0.6
    gts.to_MAC()
    expected[:, 1, :] = expected[:, 1, ::-1]
    np.testing.assert_allclose(gts.data, expected)
    assert gts.variants["maf"][1] == 0.4

    # try to do the MAC conversion again - it should fail b/c we've already done it
    with pytest.raises(AssertionError):
        gts.to_MAC()

def test_load_phenotypes():
    # create a phenotype vector with shape: samples
    expected = np.array([1, 1, 2, 2, 0])

    # can we load the data from the phenotype file?
    phens = Phenotypes(DATADIR.joinpath("simple.tsv"))
    phens.load()
    np.testing.assert_allclose(phens.data, expected)
    assert phens.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

    # try loading the data again - it should fail b/c we've already done it
    with pytest.raises(AssertionError):
        phens.load()
