from itertools import product

import pytest
import numpy as np

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
    if len(data.shape) > 1:
        data = np.squeeze(data)
    phens.data = data
    # check: are all of the phenotype values the same?
    if np.all(phens.data == phens.data[0]):
        # if they're all the same, don't standardize b/c then we'll end up with NAs
        phens.data = np.zeros(phens.data.shape)
    else:
        phens.standardize()
    return phens


def _view_tree_haps(tree) -> list:
    """
    Return the haplotype contents of a tree in an easily viewable form
    """
    return [
        [(node["label"], node["allele"]) for node in haplotype]
        for haplotype in tree.haplotypes()
    ]


def test_one_snp_perfect():
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
    tree = TreeBuilder(gens, phens).run()
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # two haplotypes: with the same SNP but different alleles
    assert len(haps) == 2
    assert len(haps[0]) == 1
    assert haps[0][0]["variant"].id == "snp0"
    assert haps[0][0]["allele"] == 0
    assert len(haps[1]) == 1
    assert haps[1][0]["variant"].id == "snp0"
    assert haps[1][0]["allele"] == 1


def test_one_snp_not_causal():
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
    phens = _create_fake_phens(np.ones(gens.data.sum(axis=2).shape))

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run()
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # no haplotypes!
    assert len(haps) == 0


def test_two_snps_single_association():
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
    phens = _create_fake_phens(gens.data.sum(axis=2)[:, 0] * 0.5)

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run()
    haps = _view_tree_haps(tree)

    # check: did the output turn out how we expected?
    # two haplotypes: one for each allele of the first SNP
    assert len(haps) == 2
    assert len(haps[0]) == 1
    assert len(haps[1]) == 1
    assert haps[0][0] == ("snp0", 0)
    assert haps[1][0] == ("snp0", 1)


@pytest.mark.xfail(reason="not implemented yet")
def test_two_snps_independent_perfect():
    """
    Two independent causal SNPs with perfect phenotype associations
    Y = 0.5 * X1 + 0.5 * X2
    """
    split_list_in_half = lambda pair: [pair[:2], pair[2:]]
    gens = _create_fake_gens(
        np.array(
            list(map(split_list_in_half, product([0, 1], repeat=4))), dtype=np.bool_
        )
    )
    gts = gens.data.sum(axis=2)
    phens = _create_fake_phens(gts[:, 0] * 0.5 + gts[:, 1] * 0.5)

    # run the treebuilder and extract the haplotypes
    builder = TreeBuilder(gens, phens, AssocTestSimple(pval_thresh=2))
    builder.run()
    tree = builder.tree
    haps = tree.haplotypes()

    # TODO: we need to handle this case, somehow
    assert False


def test_two_snps_one_branch_perfect():
    """
    Two causal SNPs on a single haplotype with perfect phenotype associations
    Y = 0.5 * ( X1 && X2 )
    This should yield two haplotypes with both SNPs having the same allele.
    The psuedocode looks like:
        if X1:
            return X2
        return 0
    """
    split_list_in_half = lambda pair: [pair[:2], pair[2:]]
    gens = np.array(
        list(map(split_list_in_half, product([0, 1], repeat=4))), dtype=np.bool_
    )
    gens = _create_fake_gens(gens)
    gts = gens.data
    phens = _create_fake_phens(0.5 * (gts[:, 0] & gts[:, 1]).sum(axis=1))

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run()
    haps = _view_tree_haps(tree)

    # check: did the output turn out how we expected?
    # one haplotype with just the minor allele of the first SNP
    # and one haplotype with both major alleles of both SNPs
    assert len(haps) == 2
    assert len(haps[0]) == 1
    assert haps[0][0] == ("snp0", 0)
    assert len(haps[1]) == 2
    assert haps[1][0] == ("snp0", 1)
    assert haps[1][1] == ("snp1", 1)


def test_two_snps_one_branch_perfect_opposite_allele():
    """
    Two causal SNPs on a single haplotype with perfect phenotype associations
    Y = 0.5 * ( X1 && X2 )
    This should yield two haplotypes with the SNPs having different alleles.
    The psuedocode looks like:
        if X1:
            return X2
        return 0
    """
    split_list_in_half = lambda pair: [pair[:2], pair[2:]]
    gens = np.array(
        list(map(split_list_in_half, product([0, 1], repeat=4))), dtype=np.bool_
    )
    gens = _create_fake_gens(gens)
    gts = gens.data
    phens = _create_fake_phens(0.5 * (gts[:, 0] & ~gts[:, 1]).sum(axis=1))

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run()
    haps = _view_tree_haps(tree)

    # check: did the output turn out how we expected?
    # one haplotype with just the minor allele of the first SNP
    # and one haplotype with a mix of alleles of both SNPs
    assert len(haps) == 2
    assert len(haps[0]) == 1
    assert haps[0][0] == ("snp0", 0)
    assert len(haps[1]) == 2
    assert haps[1][0] == ("snp0", 1)
    assert haps[1][1] == ("snp1", 0)


def test_two_snps_one_branch_perfect_opposite_direction():
    """
    This is the same as test_two_snps_one_branch_perfect except that the phenotype is
    associated in the opposite direction (ie a negative effect size) this time
    """
    split_list_in_half = lambda pair: [pair[:2], pair[2:]]
    gens = _create_fake_gens(
        np.array(
            list(map(split_list_in_half, product([0, 1], repeat=4))), dtype=np.bool_
        )
    )
    gts = gens.data
    phens = _create_fake_phens(-0.5 * (gts[:, 0] & ~gts[:, 1]).sum(axis=1))

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run()
    haps = _view_tree_haps(tree)

    # check: did the output turn out how we expected?
    # one haplotype with just the minor allele of the first SNP
    # and one haplotype with a mix of alleles of both SNPs
    assert len(haps) == 2
    assert len(haps[0]) == 1
    assert haps[0][0] == ("snp0", 0)
    assert len(haps[1]) == 2
    assert haps[1][0] == ("snp0", 1)
    assert haps[1][1] == ("snp1", 0)


def test_three_snps_one_branch_one_snp_not_causal():
    """
    Twp causal SNPs with perfect phenotype association and one SNP that isn't causal
    Y = 0.5 * ( X1 && X2 )
    This should yield a single haplotype with both SNPs having the same allele.
    The psuedocode looks like:
        if X1:
            return X2
        return 0
    """
    split_list_in_half = lambda pair: [pair[:2], pair[2:], [1, 1]]
    gens = _create_fake_gens(
        np.array(
            list(map(split_list_in_half, product([0, 1], repeat=4))), dtype=np.bool_
        )
    )
    gts = gens.data
    phens = _create_fake_phens(0.5 * (gts[:, 0] & gts[:, 1]).sum(axis=1))

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run()
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # one haplotype: with one SNP
    assert len(haps) == 1
    assert len(haps[0]) == 2
    assert haps[0][1]["variant"].id == "snp0"
    assert haps[0][2]["variant"].id == "snp1"
    assert haps[0][2]["allele"] == 1


@pytest.mark.xfail(reason="not implemented yet")
def test_four_snps_two_independent_trees_perfect():
    """
    Two independent causal SNPs each sharing a haplotype with another, different SNP
    via perfect phenotype associations
    Y = 0.5 * ( X1 && X3 ) + 0.5 * ( X2 && X4 )
    This should yield two haplotypes from different trees, where X3 occurs in the first
    and X4 occurs in the second
    """
    # a function for splitting a list into a list of pairs
    split_list = lambda pair: [pair[i : i + 2] for i in range(0, len(pair), 2)]
    gens = _create_fake_gens(
        np.array(list(map(split_list, product([0, 1], repeat=8))), dtype=np.bool_)
    )
    gts = gens.data
    phens = _create_fake_phens(
        0.5 * (gts[:, 0] & gts[:, 2]).sum(axis=1)
        + 0.5 * (gts[:, 1] & gts[:, 3]).sum(axis=1)
    )

    # TODO: we need to handle this case, somehow
    # # run the treebuilder and extract the haplotypes
    # builder = TreeBuilder(gens, phens, AssocTestSimple(pval_thresh=2))
    # builder.run()
    # tree = builder.tree
    # haps = tree.haplotypes()
    assert False


@pytest.mark.xfail(reason="not implemented yet")
def test_four_snps_two_independent_trees_perfect_one_snp_not_causal():
    """
    Two independent causal SNPs each sharing a haplotype with another, different SNP
    via perfect phenotype associations
    Y = 0.5 * ( X1 && X3 ) + 0.5 * ( X2 && X4 )
    This should yield two haplotypes from different trees, where X3 occurs in the first
    and X4 occurs in the second
    """
    # a function for splitting a list into a list of pairs
    split_list = lambda pair: [pair[i : i + 2] for i in range(0, len(pair), 2)] + [
        [1, 1]
    ]
    gens = _create_fake_gens(
        np.array(list(map(split_list, product([0, 1], repeat=8))), dtype=np.bool_)
    )
    gts = gens.data
    phens = _create_fake_phens(
        0.5 * (gts[:, 0] & gts[:, 2]).sum(axis=1)
        + 0.5 * (gts[:, 1] & gts[:, 3]).sum(axis=1)
    )

    # TODO: we need to handle this case, somehow
    # # run the treebuilder and extract the haplotypes
    # builder = TreeBuilder(gens, phens, AssocTestSimple(pval_thresh=2))
    # builder.run()
    # tree = builder.tree
    # haps = tree.haplotypes()
    assert False


@pytest.mark.xfail(reason="not implemented yet")
def test_four_snps_two_independent_trees_perfect_two_snps_not_causal():
    """
    Two independent causal SNPs each sharing a haplotype with another, different SNP
    via perfect phenotype associations plus an extra non-causal SNP (so five SNPs total)
    Y = 0.5 * ( X1 && X3 ) + 0.5 * ( X2 && X4 )
    This should yield two haplotypes from different trees, where X3 occurs in the first
    and X4 occurs in the second
    """
    # a function for splitting a list into a list of pairs
    split_list = lambda pair: [pair[i : i + 2] for i in range(0, len(pair), 2)] + [
        [1, 1],
        [0, 0],
    ]
    gens = _create_fake_gens(
        np.array(list(map(split_list, product([0, 1], repeat=8))), dtype=np.bool_)
    )
    gts = gens.data
    phens = _create_fake_phens(
        0.5 * (gts[:, 0] & gts[:, 2]).sum(axis=1)
        + 0.5 * (gts[:, 1] & gts[:, 3]).sum(axis=1)
    )

    # TODO: we need to handle this case, somehow
    # # run the treebuilder and extract the haplotypes
    # builder = TreeBuilder(gens, phens, AssocTestSimple(pval_thresh=2))
    # builder.run()
    # tree = builder.tree
    # haps = tree.haplotypes()
    assert False


@pytest.mark.xfail(reason="not implemented yet")
def test_three_snps_two_independent_trees_perfect():
    """
    Two independent causal SNPs each sharing a haplotype with the third SNP via
    perfect phenotype associations
    Y = 0.5 * ( X1 && X3 ) + 0.5 * ( X2 && X3 )
    This should yield two haplotypes from different trees, where X3 occurs in both but
    X1 occurs only in one and X2 occurs only in the other
    """
    split_list = lambda pair: [pair[:2], pair[2:4], pair[4:]]
    gens = _create_fake_gens(
        np.array(list(map(split_list, product([0, 1], repeat=6))), dtype=np.bool_)
    )
    gts = gens.data
    phens = _create_fake_phens(
        0.5 * (gts[:, 0] & gts[:, 2]).sum(axis=1)
        + 0.5 * (gts[:, 1] & gts[:, 2]).sum(axis=1)
    )

    # TODO: we need to handle this case, somehow
    # # run the treebuilder and extract the haplotypes
    # builder = TreeBuilder(gens, phens, AssocTestSimple(pval_thresh=2))
    # builder.run()
    # tree = builder.tree
    # haps = tree.haplotypes()
    assert False


@pytest.mark.xfail(reason="not implemented yet")
def test_three_snps_two_independent_trees_perfect_one_snp_not_causal():
    """
    Two independent causal SNPs each sharing a haplotype with the third SNP via
    perfect phenotype associations plus an extra non-causal SNP (so 4 SNPs total)
    Y = 0.5 * ( X1 && X3 ) + 0.5 * ( X2 && X3 )
    This should yield two haplotypes from different trees, where X3 occurs in both but
    X1 occurs only in one and X2 occurs only in the other
    """
    split_list = lambda pair: [pair[:2], pair[2:4], pair[4:]] + [[1, 1]]
    gens = _create_fake_gens(
        np.array(list(map(split_list, product([0, 1], repeat=6))), dtype=np.bool_)
    )
    gts = gens.data
    phens = _create_fake_phens(
        0.5 * (gts[:, 0] & gts[:, 2]).sum(axis=1)
        + 0.5 * (gts[:, 1] & gts[:, 2]).sum(axis=1)
    )

    # TODO: we need to handle this case, somehow
    # # run the treebuilder and extract the haplotypes
    # builder = TreeBuilder(gens, phens, AssocTestSimple(pval_thresh=2))
    # builder.run()
    # tree = builder.tree
    # haps = tree.haplotypes()
    assert False


def test_two_snps_two_branches_perfect():
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
    gens = _create_fake_gens(
        np.array(
            list(map(split_list_in_half, product([0, 1], repeat=4))), dtype=np.bool_
        )
    )
    gts = gens.data
    phens = _create_fake_phens(0.5 * (gts[:, 0] | gts[:, 1]).sum(axis=1))

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run()
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # two haplotypes, each with both SNPs
    assert len(haps) == 2
    assert len(haps[0]) == 2
    for i in range(2):
        assert haps[i][1]["variant"].id == "snp0"
        assert haps[i][2]["variant"].id == "snp1"
        assert haps[i][2]["allele"] == 1


def test_two_snps_two_branches_perfect_one_snp_not_causal():
    """
    Two causal SNPs on a single haplotype with perfect phenotype associations
    plus an extra non-causal SNP
    Y = 0.5 * ( X1 || X2 )
    This should yield two different haplotypes for each of the alleles.
    The pseudocode looks like:
        if X1:
            return 1
        else:
            return X2
    """
    split_list_in_half = lambda pair: [pair[:2], pair[2:]] + [[1, 1]]
    gens = _create_fake_gens(
        np.array(
            list(map(split_list_in_half, product([0, 1], repeat=4))), dtype=np.bool_
        )
    )
    gts = gens.data
    phens = _create_fake_phens(0.5 * (gts[:, 0] | gts[:, 1]).sum(axis=1))

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens).run()
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # two haplotypes, each with both SNPs
    assert len(haps) == 2
    assert len(haps[0]) == 2
    for i in range(2):
        assert haps[i][1]["variant"].id == "snp0"
        assert haps[i][2]["variant"].id == "snp1"
        assert haps[i][2]["allele"] == 1


def test_ppt_case():
    """
    Test the example case from my powerpoint slides
    This is a more complicated example.
    Two causal SNPs on a single haplotype and three causal SNPs on another haplotype,
    both with perfect phenotype associations plus an extra non-causal SNP
    The two haplotypes occur on the same tree and share two SNPs: one being the root
    and another being SNP2 (with different alleles)
    Y = 0.5 * ( ( X1 && X2 ) || ( X1 && X2 && X3 ) )
    This should yield two different haplotypes for each of the alleles.
    The pseudocode looks like:
        if X1:
            return X2
        else:
            if X2:
                return 0
            return X3
    """
    # a function for splitting a list into a list of pairs
    split_list = lambda pair: [pair[i : i + 2] for i in range(0, len(pair), 2)]
    # create genotypes for 3 samples, 4 SNPs
    gens = _create_fake_gens(
        np.array(
            list(map(split_list_in_half, product([0, 1], repeat=4 * 2))), dtype=np.bool_
        )
    )
    gts = gens.data
    # create phenotypes for 3 samples, excluding the last
    phens = _create_fake_phens(
        0.5
        * ((gts[:, 0] & gts[:, 1]) | (~gts[:, 0] & ~gts[:, 1] & gts[:, 2])).sum(axis=1)
    )

    # run the treebuilder and extract the haplotypes
    tree = TreeBuilder(gens, phens, AssocTestSimple(pval_thresh=2)).run()
    haps = tree.haplotypes()

    # check: did the output turn out how we expected?
    # two haplotypes: one with three SNPs and one with two
    assert len(haps) == 2
    assert tuple([len(hap) for hap in haps]) == (3, 2)
    for i in range(3):
        assert haps[0][i]["variant"].id == "snp" + str(i)
    for i in range(2):
        assert haps[1][i]["variant"].id == "snp" + str(i)
