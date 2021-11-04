from happler.data import Data


def test_load_genotypes():
	assert Data('tests/example.vcf', '') == np.zeros(20).reshape((4,5))
