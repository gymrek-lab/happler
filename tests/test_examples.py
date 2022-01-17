import numpy as np

from itertools import product

from happler.data import Genotypes, Phenotypes
from happler.tree import TreeBuilder, AssocTestSimple


def _create_fake_gens(data) -> Genotypes:
    """
    Create a fake Genotypes object for testing purposes

    Parameters
    ----------
    data : np.array
        The genotype matrix with shape n x p x 2 of dtype np.bool_

    Returns
    -------
    Genotypes
        The genotypes object with sample names and variant info filled in
    """
    gens = Genotypes(fname=None)
    gens.samples = tuple("samp" + str(i) for i in range(data.shape[0]))
    gens.variants = np.array(
        [("snp" + str(i), "chr0", i, 0.75) for i in range(data.shape[1])],
        dtype=[
            ("id", "U50"),
            ("chrom", "U10"),
            ("pos", np.uint),
            ("aaf", np.float64),
        ],
    )
    gens.data = data
    return gens


def _create_fake_phens(data) -> Phenotypes:
    """
    Create a fake Phenotypes object for testing purposes

    Parameters
    ----------
    data : np.array
        The phenotypes matrix with shape n x 1 of dtype np.float64

    Returns
    -------
    Phenotypes
        The phenotypes object with sample names filled in
    """
    phens = Phenotypes(fname=None)
    phens.samples = tuple("samp" + str(i) for i in range(data.shape[0]))
    phens.data = data
    phens.standardize()
    return phens


def test_treebuilder_one_snp_perfect():
    """
    The most simple case. One causal SNP with a perfect phenotype association:
    Y = 0.5 * X1
    """
    # 4 samples
    gens = _create_fake_gens(
        np.array(
            [
                [[0, 0]],
                [[0, 1]],
                [[1, 0]],
                [[1, 1]],
            ],
            dtype=np.bool_,
        )
    )
    phens = _create_fake_phens(gens.data.sum(axis=2) * 0.5)

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run(root=0)
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # one haplotype: with one SNP
    assert len(haps) == 1
    assert len(haps[0]) == 1
    assert haps[0][0]["variant"].id == "snp0"


def test_treebuilder_one_snp_not_causal():
    """
    One non-causal SNP with no phenotype association
    """
    # 4 samples
    gens = _create_fake_gens(
        np.array(
            [
                [[0, 0]],
                [[0, 1]],
                [[1, 0]],
                [[1, 1]],
            ],
            dtype=np.bool_,
        )
    )
    phens = _create_fake_phens(np.ones(gens.data.sum(axis=2).shape) * 0.5)

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run(root=0)
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # one haplotype: with one SNP
    assert len(haps) == 1
    assert len(haps[0]) == 1
    assert haps[0][0]["variant"].id == "snp0"


def test_treebuilder_two_snps_single_association():
    """
    One causal SNP with a perfect phenotype association and one SNP that isn't causal
    Y = 0.5 * X1
    """
    # 4 samples
    gens = _create_fake_gens(
        np.array(
            [
                [[0, 0], [1, 1]],
                [[0, 1], [1, 1]],
                [[1, 0], [1, 1]],
                [[1, 1], [1, 1]],
            ],
            dtype=np.bool_,
        )
    )
    phens = _create_fake_phens(gens.data[:,0].sum(axis=2) * 0.5)

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run(root=0)
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # one haplotype: with one SNP
    assert len(haps) == 1
    assert len(haps[0]) == 1
    assert haps[0][0]["variant"].id == "snp0"


def test_treebuilder_two_snps_independent_perfect():
    """
    Two independent causal SNPs with perfect phenotype associations
    Y = 0.5 * X1 + 0.5 * X2
    """
    split_list_in_half = lambda pair: [pair[:2], pair[2:]]
    # 4 samples
    gens = _create_fake_gens(
        np.array(
            list(map(split_list_in_half, product([0, 1], repeat=4))), dtype=np.bool_
        )
    )
    gts = gens.data.sum(axis=2)
    phens = _create_fake_phens(gts[:, 0] * 0.5 + gts[:, 1] * 0.5)

    # TODO: we need to handle this case, somehow
    # # run the treebuilder and extract the haplotypes
    # builder = TreeBuilder(gens, phens, AssocTestSimple(pval_thresh=2))
    # builder.run(root=0)
    # tree = builder.tree
    # haps = tree.haplotypes()


def test_treebuilder_two_snps_one_branch_perfect():
    """
    Two causal SNPs on a single haplotype with perfect phenotype associations
    Y = 0.5 * ( X1 && X2 )
    This should yield a single haplotype with both SNPs having the same allele.
    The psuedocode looks like:
        if X1:
            return X2
    """
    split_list_in_half = lambda pair: [pair[:2], pair[2:]]
    # 4 samples
    gens = _create_fake_gens(
        np.array(
            list(map(split_list_in_half, product([0, 1], repeat=4))), dtype=np.bool_
        )
    )
    gts = gens.data
    phens = _create_fake_phens(0.5 * (gts[:, 0] & gts[:, 1]).sum(axis=1))

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run(root=0)
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # one haplotype with both SNPs
    assert len(haps) == 1
    assert len(haps[0]) == 2
    assert haps[0][1]["variant"].id == "snp0"
    assert haps[0][2]["variant"].id == "snp1"
    assert haps[0][2]["allele"] == 1


def test_treebuilder_two_snps_two_branches_perfect():
    """
    Two causal SNPs on a single haplotype with perfect phenotype associations
    Y = 0.5 * ( X1 || X2 )
    This should yield two different haplotypes for each of the alleles.
    The pseudocode looks like:
        if X1:
            return 1
        else:
            return X2

    """
    split_list_in_half = lambda pair: [pair[:2], pair[2:]]
    # 4 samples
    gens = _create_fake_gens(
        np.array(
            list(map(split_list_in_half, product([0, 1], repeat=4))), dtype=np.bool_
        )
    )
    gts = gens.data
    phens = _create_fake_phens(0.5 * (gts[:, 0] | gts[:, 1]).sum(axis=1))

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run(root=0)
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # two haplotypes, each with both SNPs
    assert len(haps) == 2
    assert len(haps[0]) == 2
    for i in range(2):
        assert haps[i][1]["variant"].id == "snp0"
        assert haps[i][2]["variant"].id == "snp1"
        assert haps[i][2]["allele"] == 1


def test_treebuilder_ppt_case():
    """
    Test the example case from my powerpoint slides
    This is a more complicated example.
    TODO: adjust the genotypes and phenotypes to achieve the phenotypes I want
    """
    # create genotypes for 3 samples, 4 SNPs
    gens = _create_fake_gens(
        np.array(
            [
                [[0, 0], [1, 1], [1, 1], [1, 1]],
                [[0, 1], [1, 0], [0, 0], [1, 0]],
                [[1, 1], [0, 0], [0, 0], [0, 0]],
            ],
            dtype=np.bool_,
        )
    )
    # create phenotypes for 3 samples
    phens = _create_fake_phens(np.array([3, 6, 7], dtype=np.float64))

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens, AssocTestSimple(pval_thresh=2)).run(root=0)
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # two haplotypes: one with three SNPs and one with two
    assert len(haps) == 2
    assert tuple([len(hap) for hap in haps]) == (3, 2)
    for i in range(3):
        assert haps[0][i]["variant"].id == "snp" + str(i)
    for i in range(2):
        assert haps[1][i]["variant"].id == "snp" + str(i)
