from __future__ import annotations
from logging import Logger

import numpy as np
from haptools.logging import getLogger
from haptools.data import Genotypes, Phenotypes

from .variant import Variant
from .haplotypes import Haplotype
from .tree import Tree, NodeResults
from .assoc_test import AssocTest, AssocTestSimple
from .terminator import Terminator, TTestTerminator, BICTerminator


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
    method: AssocTest, optional
        The type of association test to perform at each node when constructing the tree
    terminator: Terminator, optional
        The type of test to use for deciding whether to terminate a branch
    results: type[NodeResults]
        The class to use when instantiating the results of a node association test
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    >>> gens = Genotypes.load('tests/data/simple.vcf')
    >>> phens = Phenotypes.load('tests/data/simple.tsv')
    >>> tree = TreeBuilder(gens, phens).run()
    """

    def __init__(
        self,
        genotypes: Genotypes,
        phenotypes: Phenotypes,
        method: AssocTest = AssocTestSimple(),
        terminator: Terminator = TTestTerminator(),
        results_type: type[NodeResults] = NodeResults,
        log: Logger = None,
    ):
        self.gens = genotypes
        self.phens = phenotypes
        self.method = method
        self.terminator = terminator
        self.results_type = results_type
        self.tree = None
        self.log = log or getLogger(self.__class__.__name__)

    def __repr__(self):
        return str((self.gens, self.phens))

    def run(self):
        """
        Run the tree builder and create a tree rooted at the provided variant
        """
        if self.tree is not None:
            raise AssertionError(
                "A tree already exists for this TreeBuilder. Please create a new one."
            )
        # step one: initialize the tree
        self.tree = Tree()
        # step two: create a haplotype
        # we're at the root of tree, so we need to create a new, empty haplotype
        parent_hap = Haplotype(num_samples=len(self.gens.samples))
        # step three: create the rest of the tree
        self._create_tree(parent_hap, parent_idx=0)
        return self.tree

    def _create_tree(
        self,
        parent_hap: Haplotype,
        parent_idx: int,
        parent_res: NodeResults = None,
    ):
        """
        Recursive helper to the run() function

        Test for an association with the parent_hap and its allele, and then add a
        subtree under the node with index parent_idx in the tree

        Parameters
        ----------
        parent_hap : Haplotype
            The haplotype containing all variants up to (but NOT including) the parent
        parent_idx : int
            The index of the parent node in the tree
        parent_res : NodeResults
            The results of the association test for the parent node
        """
        if len(parent_hap.nodes):
            parent = parent_hap.nodes[-1]
            # add the best variant as a node in the tree
        self.log.debug(
            "Adding variants to "
            + (
                "parent {} with allele {}".format(parent[0].id, parent[1])
                if len(parent_hap.nodes)
                else "root"
            )
        )
        # find the variant-allele pairs that give the best haplotype
        for variant, allele, results in self._find_split(parent_hap, parent_res):
            if variant is None:
                # there were no significant variants!
                continue
            new_node_idx = self.tree.add_node(variant, parent_idx, allele, results)
            # create a new Haplotype with the variant-allele pair added
            variant_gts = self.gens.data[:, variant.idx, :2] == allele
            new_parent_hap = parent_hap.append(variant, allele, variant_gts)
            self._create_tree(new_parent_hap, new_node_idx, results)

    def _find_split(
        self, parent: Haplotype, parent_res: NodeResults = None
    ) -> tuple[Variant, np.void]:
        """
        Find the variant/allele that best fits under the parent_idx node

        Parameters
        ----------
        parent : Haplotype
            The haplotype containing all variants up to (and including) the parent
        parent_res : NodeResults, optional
            The results of the tests performed on the parent node

        Returns
        -------
        tuple[Variant, float]
            The variant that best fits under the parent node with the allele edge AND
            the results (ex: beta, pval) of the haplotype association test after
            incorporating that variant
        """
        num_samps = len(self.gens.samples)
        results = {}
        best_p_idx = {}
        # iterate through the two possible alleles and try all SNPs with that allele
        alleles = (0, 1)
        for allele in alleles:
            # step 1: transform the GT matrix into a haplotype matrix
            hap_matrix = parent.transform(self.gens, allele)
            if hap_matrix.shape[1] == 0:
                # if there weren't any genotypes left, just return None
                yield None, allele, None
                continue
            # step 2: run all association tests on all of the haplotypes
            results[allele] = self.method.run(
                hap_matrix.sum(axis=2),
                self.phens.data[:, 0],
            )
            # also, record the best p-value among all the SNPs with this allele
            best_p_idx[allele] = results[allele].data["pval"].argmin()
        # exit if neither of the alleles worked
        if not len(best_p_idx):
            return
        # step 3: find the index of the best variant within the haplotype matrix
        best_allele = min(
            best_p_idx, key=lambda a: results[a].data["pval"][best_p_idx[a]]
        )
        best_var_idx = best_p_idx[best_allele]
        num_tests = len(parent.nodes) + 1
        # step 4: find the index of the best variant within the genotype matrix
        # we need to account for indices that we removed when running transform()
        # There might be a faster way of doing this but for now we're just going to
        # live with it
        for gt_idx in sorted(parent.node_indices):
            # add each idx back in, so long as they are less than the target idx
            if gt_idx > best_var_idx:
                break
            best_var_idx += 1
        # step 5: retrieve the Variant with the best p-value
        best_variant = Variant.from_np(self.gens.variants[best_var_idx], best_var_idx)
        self.log.debug("Chose variant {}".format(best_variant.id))
        # iterate through all of the alleles of the best variant and check if they're
        # significant
        for allele in alleles:
            best_results = results[allele].data[best_p_idx[best_allele]]
            node_res = self.results_type.from_np(best_results)
            # step 6: check whether we should terminate the branch
            self.log.debug(
                "Testing variant {} / allele {} with parent_res {} and node_res {}"
                .format(best_variant.id, allele, parent_res, node_res)
            )
            if self.terminator.check(parent_res, node_res, num_samps, num_tests):
                yield None, allele, node_res
                continue
            yield best_variant, allele, node_res
