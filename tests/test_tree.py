import pytest
import numpy as np
from pathlib import Path

from happler.data import VariantType, Genotypes, Phenotypes
from happler.tree import Node, Tree, TreeBuilder


DATADIR = Path(__file__).parent.joinpath("data")


def test_node():
    node = Node(idx=0, ID="SNP0", pos=1)
    assert node.type == VariantType("snp")
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
    assert node == Node.from_np(node.idx, np_node[0])


def test_tree():
    node = Node(idx=0, ID="SNP0", pos=1)
    tree = Tree(root_node=node)


def test_tree_builder():
    gens = Genotypes.load(DATADIR.joinpath("simple.vcf.gz"))
    phens = Phenotypes.load(DATADIR.joinpath("simple.tsv"))
    tree_builder = TreeBuilder(gens, phens)
