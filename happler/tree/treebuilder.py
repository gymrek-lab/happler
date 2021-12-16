from .tree import Tree
from ..data import Genotypes, Phenotypes
from .haplotypes import Variant, Haplotype
from .test_assoc import TestAssoc, TestAssocSimple


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
    method: TestAssoc
        The type of association test to perform at each node when constructing the tree
    """

    def __init__(self, genotypes: Genotypes, phenotypes: Phenotypes, method: TestAssoc):
        self.gens = genotypes
        self.phens = phenotypes
        self.method = method
        self.tree = None

    def __repr__(self):
        return self.gens, self.phens

    def run(self, root: int):
        """
        Run the tree builder and create a tree rooted at the provided variant

        Parameters
        ----------
        root : int
            The index of the variant to use at the root of tree. This should be an
            index into :py:attr:`~.TreeBuilder.gens`
        """
        # step one: initialize the tree
        root_node = Variant.from_np(self.gens[root], idx=root)
        self.tree = Tree(root_node)
        # step two: create the tree
        self._create_tree(root_node)

    def _create_tree(
        self, parent: Variant, parent_hap: Haplotype = None, parent_idx: int = 0
    ):
        """
        Recursive helper to the run() function

        Adds a subtree under the node with index parent_idx in the tree

        Parameters
        ----------
        parent : Variant
            An existing node in the tree under which we should consider creating a subtree
        parent_hap : Haplotype
            The haplotype containing all variants up to (but NOT including) the parent
        parent_idx : int
            The index of the parent node in the tree
        """
        # we consider two possible alleles
        alleles = (0, 1)
        for allele in alleles:
            # find the best variant, add it to the tree, and then create a new subtree
            # under it
            best_variant = self._find_split(parent_idx, allele)
            new_node_idx = self.tree.add_node(best_variant, parent_idx, allele)
            if parent_hap is None:
                parent_hap = Haplotype()
            parent_hap.append(parent, allele)
            self._create_tree(best_variant, parent_hap, new_node_idx)

    def _find_split(self, parent: Haplotype, parent_idx: int, allele: int) -> Variant:
        """
        Find the variant that best fits under the parent_idx node with the allele edge

        Parameters
        ----------
        TODO
        """
        # step 1: transform GT matrix into haplotype matrix
        hap_matrix = parent.transform(self.gens)
        # step 2: test assoc
        p_values = self.method.run(hap_matrix, self.phenotypes.data)
        # step 3:
