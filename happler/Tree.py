import networkx as nx

from .Data import Genotypes, Phenotypes


class Tree:
    """
    A tree where

    1. nodes are variants (and also haplotypes)
    2. each node can have at most two branches for each of its alleles

    Attributes
    ----------
    graph: nx.DiGraph
        The underlying directed acyclic graph representing the tree
    """

    def __init__(self):
        self.graph = nx.DiGraph()


class TreeBuilder:
    """
    Creates a Tree object from provided Genotypes and Phenotypes

    Attributes
    ----------
    tree: Tree
        A tree representing haplotypes composed from trait-associated genotypes
    gens: Genotypes
        The genotypes from which the tree should be built
    phens: Phenotypes
        The phenotypes from which the tree should be built
    """

    def __init__(self, genotypes: Genotypes, phenotypes: Phenotypes):
        self.gens = genotypes
        self.phens = phenotypes
        self.tree = Tree()

    def __repr__(self):
        return self.gens, self.phens
