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
    NodeResults,
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
    tree = Tree()

    snp0 = Variant(idx=0, id="SNP0", pos=1)
    tree.add_node(snp0, parent_idx=0, allele=0)
    assert tree.num_nodes, tree.num_variants == (2, 1)
    tree.add_node(snp0, parent_idx=0, allele=1)
    assert tree.num_nodes, tree.num_variants == (3, 1)

    snp1 = Variant(idx=1, id="SNP1", pos=2)
    with pytest.raises(ValueError):
        tree.add_node(snp1, parent_idx=0, allele=0)
    assert tree.num_nodes, tree.num_variants == (3, 1)
    tree.add_node(snp1, parent_idx=1, allele=0)
    assert tree.num_nodes, tree.num_variants == (4, 2)


def test_tree_dot():
    node = Variant(idx=0, id="SNP0", pos=1)
    tree = Tree()
    node_idx = tree.add_node(node, 0, 0)
    node_idx = tree.add_node(Variant(idx=1, id="SNP1", pos=2), node_idx, 0)
    node_idx = tree.add_node(Variant(idx=2, id="SNP2", pos=3), node_idx, 1)
    snp3 = Variant(idx=3, id="SNP3", pos=4)
    tree.add_node(snp3, node_idx, 0)
    tree.add_node(snp3, node_idx, 1)
    assert (
        tree.dot()
        == "strict digraph  {\n0 [label=root];\n1 [label=SNP0];\n2 [label=SNP1];\n3"
        " [label=SNP2];\n4 [label=SNP3];\n5 [label=SNP3];\n0 -> 1  [label=0];\n1 ->"
        " 2  [label=0];\n2 -> 3  [label=1];\n3 -> 4  [label=0];\n3 -> 5 "
        " [label=1];\n}\n"
    )


def test_get_haplotypes_from_tree():
    # create three new nodes
    snp1 = Variant(idx=0, id="SNP1", pos=1)
    snp2 = Variant(idx=1, id="SNP2", pos=2)
    snp3 = Variant(idx=2, id="SNP3", pos=3)

    # create a tree composed of just one variant
    tree = Tree()
    snp1_idx = tree.add_node(snp1, parent_idx=0, allele=0)
    haps = tree.haplotypes()
    assert len(haps) == 1
    assert haps[0][0]["variant"] == snp1

    # add two child nodes
    snp2_idx = tree.add_node(snp2, parent_idx=snp1_idx, allele=0)
    tree.add_node(snp2, parent_idx=snp1_idx, allele=1)
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


@pytest.mark.skip(reason="test not completely written yet")
def test_tree_builder():
    gens = Genotypes.load(DATADIR.joinpath("simple.vcf"))
    phens = Phenotypes.load(DATADIR.joinpath("simple.tsv"))
    tree_builder = TreeBuilder(gens, phens)
    # tree = tree_builder.run()
    # TODO: assert that the tree looks correct
    assert False


def test_haplotype():
    node = Variant(idx=0, id="SNP0", pos=1)
    gts = np.array(
        [
            [[0, 0], [1, 1], [1, 0]],
            [[0, 1], [1, 0], [0, 0]],
            [[1, 0], [0, 0], [1, 1]],
            [[1, 1], [0, 1], [0, 0]],
        ],
        dtype=np.bool_,
    )
    node_allele = 0
    hap_data = gts[:, node.idx, :] == node_allele
    hap = Haplotype(((node, node_allele),), hap_data)
    assert hap.nodes == Haplotype.from_node(node, 0, hap.data).nodes

    new_node = Variant(idx=1, id="SNP1", pos=1)
    new_node_allele = 1
    new_hap_data = gts[:, new_node.idx, :] == new_node_allele
    hap = hap.append(new_node, new_node_allele, new_hap_data)
    nodes = ((node, node_allele), (new_node, new_node_allele))
    assert hap.nodes == Haplotype(nodes, hap.data).nodes
    np.testing.assert_allclose(hap.data, hap_data & new_hap_data)


def test_haplotype_transform():
    gens = Genotypes.load(DATADIR.joinpath("simple.vcf"))
    variant_idx = 0
    allele = 0
    gens.data[:, variant_idx, allele] = 1
    variant_gts = gens.data[:, variant_idx, :] == allele
    variant = Variant.from_np(gens.variants[variant_idx], variant_idx)
    hap = Haplotype.from_node(variant, allele, variant_gts)
    gens_without_variant = (
        gens.data[:, (variant_idx + 1) :, :] & variant_gts[:, np.newaxis]
    )
    np.testing.assert_allclose(hap.transform(gens), gens_without_variant)


def test_haplotypes_write():
    # create three new nodes
    snp1 = Variant(idx=0, id="SNP1", pos=1)
    snp2 = Variant(idx=1, id="SNP2", pos=2)
    snp3 = Variant(idx=2, id="SNP3", pos=3)

    # create a results object that all of the SNPs can share
    res = NodeResults(beta=0.1, pval=0.1, stderr=0.1)

    # create a tree composed of these nodes
    tree = Tree()
    snp1_idx = tree.add_node(snp1, parent_idx=0, allele=0, results=res)
    snp2_idx = tree.add_node(snp2, parent_idx=snp1_idx, allele=1, results=res)
    tree.add_node(snp3, parent_idx=snp2_idx, allele=1, results=res)
    snp1_idx = tree.add_node(snp1, parent_idx=0, allele=1, results=res)
    snp2_idx = tree.add_node(snp2, parent_idx=snp1_idx, allele=0, results=res)

    # write the tree to a file
    with open("test_write.haps", "w") as file:
        Haplotypes.from_tree(tree).write(file)

    # verify that the results appear as intended
    with open("test_write.haps", "r") as file:
        lines = [line.rstrip() for line in file.readlines()]
        assert lines == [
            "H\t0\t0\t0.00\t0.00\tnan",
            "V\tSNP1\t0\t0.10",
            "V\tSNP2\t1\t0.10",
            "V\tSNP3\t1\t0.10",
            "H\t1\t0\t0.00\t0.00\tnan",
            "V\tSNP1\t1\t0.10",
            "V\tSNP2\t0\t0.10",
        ]

    # remove the file
    os.remove("test_write.haps")


def test_simple_assoc():
    # TODO
    pass
