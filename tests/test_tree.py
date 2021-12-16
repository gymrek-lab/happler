import pytest
import numpy as np

from pathlib import Path

from happler.data import Genotypes, Phenotypes
from happler.tree import VariantType, Variant, Haplotype, Haplotypes, Tree, TreeBuilder


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

    tree.add_node(Variant(idx=0, id="SNP1", pos=2), 0, 0)

    assert tree.num_nodes == 2


def test_tree_builder():
    gens = Genotypes.load(DATADIR.joinpath("simple.vcf"))
    phens = Phenotypes.load(DATADIR.joinpath("simple.tsv"))
    tree_builder = TreeBuilder(gens, phens)
    # root_node =
    # tree_builder.run(root_node)


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


def test_haplotypes():
    gens = Genotypes.load(DATADIR.joinpath("simple.vcf"))
    haps = Haplotypes(gens)

    variant = Variant.from_np(gens.variants[0], 0)
    haps.add_variant(variant, 0)
    np.testing.assert_allclose(gens.data[:, 0, 0][:, np.newaxis], haps.data)

    variant = Variant.from_np(gens.variants[1], 1)
    haps.add_variant(variant, 0)
    np.testing.assert_allclose(
        np.logical_and(gens.data[:, 0, 0], gens.data[:, 1, 0])[:, np.newaxis], haps.data
    )

def test_simple_assoc():
    pass
