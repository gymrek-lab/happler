import pytest
import numpy as np
from pathlib import Path

from happler.data import Genotypes, Phenotypes
from happler.tree import Tree, TreeBuilder


DATADIR = Path(__file__).parent.joinpath("data")


def test_tree():
    tree = Tree()


def test_tree_builder():
    gens = Genotypes.load(DATADIR.joinpath("simple.vcf.gz"))
    phens = Phenotypes.load(DATADIR.joinpath("simple.tsv"))
    tree_builder = TreeBuilder(gens, phens)
