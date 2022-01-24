import os
from pathlib import Path

import pytest
import numpy as np

from happler.data import Genotypes, Phenotypes
from happler.tree import (
    VariantType,
    Variant,
    Haplotype,
    Haplotypes,
    Tree,
    TreeBuilder,
    AssocTestSimple,
)


DATADIR = Path(__file__).parent.joinpath("data")


def test_variant_type():
    assert VariantType("snp") == VariantType("SNP")

    assert str(VariantType("snp")) == "SNP"

    with pytest.raises(ValueError):
        VariantType("indel")


def test_node():
    node = Variant(idx=0, id="SNP0", pos=1)
    assert node.id == "SNP0"

    np_node = np.array(
        [(node.id, "1", node.pos, "0.25")],
        dtype=[
            ("id", "U50"),
            ("chrom", "U10"),
            ("pos", np.uint),
            ("aaf", np.float64),
        ],
    )
    assert node == Variant.from_np(np_node[0], node.idx)


def test_tree():
    node = Variant(idx=0, id="SNP0", pos=1)
    tree = Tree(root=node)

    assert tree.num_nodes == 1

    tree.add_node(Variant(idx=1, id="SNP1", pos=2), 0, 0)

    assert tree.num_nodes == 2


def test_tree_dot():
    node = Variant(idx=0, id="SNP0", pos=1)
    tree = Tree(root=node)
    node_idx = tree.add_node(Variant(idx=1, id="SNP1", pos=2), 0, 0)
    tree.add_node(Variant(idx=2, id="SNP2", pos=3), 0, 1)
    tree.add_node(Variant(idx=3, id="SNP3", pos=4), node_idx, 0)
    tree.add_node(Variant(idx=4, id="SNP4", pos=5), node_idx, 1)
    assert (
        tree.dot()
        == "strict digraph  {\n0 [label=SNP0];\n1 [label=SNP1];\n2 [label=SNP2];\n3"
        " [label=SNP3];\n4 [label=SNP4];\n0 -> 1  [label=0];\n0 -> 2  [label=1];\n1"
        " -> 3  [label=0];\n1 -> 4  [label=1];\n}\n"
    )


def test_get_haplotypes_from_tree():
    # create three new nodes
    snp1 = Variant(idx=0, id="SNP1", pos=1)
    snp2 = Variant(idx=1, id="SNP2", pos=2)
    snp3 = Variant(idx=2, id="SNP3", pos=3)

    # create a tree composed of just one node
    tree = Tree(root=snp1)
    haps = tree.haplotypes()
    assert len(haps) == 1
    assert haps[0][0]["variant"] == snp1

    # add two child nodes
    snp2_idx = tree.add_node(snp2, parent_idx=0, allele=0)
    tree.add_node(snp2, parent_idx=0, allele=1)
    haps = tree.haplotypes()
    assert len(haps) == 2
    assert tuple([len(hap) for hap in haps]) == (2, 2)
    assert [haps[0][i]["variant"] for i in range(2)] == [snp1, snp2]
    assert [haps[1][i]["variant"] for i in range(2)] == [snp1, snp2]

    # add a single child to the first leaf
    tree.add_node(snp3, parent_idx=snp2_idx, allele=1)
    haps = tree.haplotypes()
    assert len(haps) == 2
    assert tuple([len(hap) for hap in haps]) == (3, 2)
    assert [haps[0][i]["variant"] for i in range(3)] == [snp1, snp2, snp3]
    assert [haps[1][i]["variant"] for i in range(2)] == [snp1, snp2]


def test_tree_builder():
    gens = Genotypes.load(DATADIR.joinpath("simple.vcf"))
    phens = Phenotypes.load(DATADIR.joinpath("simple.tsv"))
    tree_builder = TreeBuilder(gens, phens)
    root_node = 0
    tree_builder.run(root_node)
    # TODO: assert that the tree looks correct


def test_haplotype():
    node = Variant(idx=0, id="SNP0", pos=1)
    gts = np.array([0, 1, 1], dtype=np.bool_)
    hap = Haplotype(((node, 0),), gts)
    assert hap.nodes == Haplotype.from_node(node, 0, hap.data).nodes

    new_node = Variant(idx=1, id="SNP1", pos=1)
    new_gts = np.array([1, 0, 1], dtype=np.bool_)
    hap = hap.append(new_node, 0, new_gts)
    nodes = ((node, 0), (new_node, 0))
    assert hap.nodes == Haplotype(nodes, hap.data).nodes
    np.testing.assert_allclose(hap.data, np.array([0, 0, 1], np.bool_))


def test_haplotype_transform():
    gens = Genotypes.load(DATADIR.joinpath("simple.vcf"))
    variant_idx = 0
    allele_idx = 0
    gens.data[:, variant_idx, allele_idx] = 1
    variant_gts = gens.data[:, variant_idx, allele_idx]
    variant = Variant.from_np(gens.variants[variant_idx], variant_idx)
    hap = Haplotype.from_node(variant, allele_idx, variant_gts)
    gens_without_variant = gens.data[:, (variant_idx + 1) :, :]
    np.testing.assert_allclose(hap.transform(gens), gens_without_variant)


def test_haplotypes_write():
    # create three new nodes
    snp1 = Variant(idx=0, id="SNP1", pos=1)
    snp2 = Variant(idx=1, id="SNP2", pos=2)
    snp3 = Variant(idx=2, id="SNP3", pos=3)

    # create a results object that all of the SNPs can share
    res = np.array([(0.1, 0.1)], dtype=[("beta", np.float64), ("pval", np.float64)])[0]

    # create a tree composed of these nodes
    tree = Tree(root=snp1)
    snp2_idx = tree.add_node(snp2, parent_idx=0, allele=0, results=res)
    tree.add_node(snp2, parent_idx=0, allele=1, results=res)
    tree.add_node(snp3, parent_idx=snp2_idx, allele=1, results=res)

    # write the tree to a file
    Haplotypes.from_tree(tree).write("test_write.haps")

    # verify that the results appear as intended
    with open("test_write.haps", "r") as file:
        lines = [line.rstrip() for line in file.readlines()]
        assert lines == [
            "H\t0\t0\t0.00\t0.00\tnan",
            "V\tSNP1\tnan\tnan",
            "V\tSNP2\t0\t0.10",
            "V\tSNP3\t1\t0.10",
            "H\t1\t0\t0.00\t0.00\tnan",
            "V\tSNP1\tnan\tnan",
            "V\tSNP2\t1\t0.10",
        ]

    # remove the file
    os.remove("test_write.haps")


def test_simple_assoc():
    # TODO
    pass
