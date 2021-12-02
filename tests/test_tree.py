import pytest
import numpy as np

from pathlib import Path

from happler.data import VariantType, Genotypes, Phenotypes
from happler.tree import Node, Tree, TreeBuilder


DATADIR = Path(__file__).parent.joinpath("data")


def test_node():
    node = Node(idx=0, id="SNP0", pos=1)
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
    assert node == Node.from_np(np_node[0], node.idx)


def test_tree():
    node = Node(idx=0, id="SNP0", pos=1)
    tree = Tree(root=node)

    assert tree.num_nodes == 1

    tree.add_node(Node(idx=0, id="SNP1", pos=2), 0, 0)

    assert tree.num_nodes == 2


def test_tree_builder():
    gens = Genotypes.load(DATADIR.joinpath("simple.vcf.gz"))
    phens = Phenotypes.load(DATADIR.joinpath("simple.tsv"))
    tree_builder = TreeBuilder(gens, phens)
    # root_node =
    # tree_builder.run(root_node)
