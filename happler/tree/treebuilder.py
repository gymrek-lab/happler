from .tree import Tree
from ..data import Genotypes, Phenotypes


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
