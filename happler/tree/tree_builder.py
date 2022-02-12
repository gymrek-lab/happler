from __future__ import annotations
from logging import getLogger, Logger

import numpy as np
from scipy.stats import t as t_dist

from .tree import Tree, NodeResults
from ..data import Genotypes, Phenotypes
from .variant import Variant
from .haplotypes import Haplotype
from .assoc_test import AssocTest, AssocTestSimple


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
        log: Logger = None,
    ):
        self.gens = genotypes
        self.phens = phenotypes
        self.method = method
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
            "Adding variants to " + (
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
            variant_gts = self.gens.data[:, variant.idx, :] == allele
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
            results[allele] = self.method.run(hap_matrix.sum(axis=2), self.phens.data)
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
        num_snps_tested = len(results[best_allele].data)
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
        self.log.debug(
            "Chose variant {}".format(best_variant.id)
        )
        # iterate through all of the alleles of the best variant and check if they're
        # significant
        for allele in alleles:
            best_results = results[allele].data[best_p_idx[best_allele]]
            node_res = NodeResults.from_np(best_results)
            # step 6: check whether we should terminate the branch
            self.log.debug(
                "Testing variant {} / allele {} with parent_res {} and node_res {}".format(
                    best_variant.id, allele, parent_res, node_res
                )
            )
            if self.check_terminate(parent_res, node_res, num_samps, num_snps_tested):
                yield None, allele, node_res
                continue
            yield best_variant, allele, node_res

    def check_terminate(
        self,
        parent_res: NodeResults,
        node_res: NodeResults,
        num_samps: int,
        num_tests: int,
    ) -> bool:
        """
        Check whether this branch should be terminated.

        Parameters
        ----------
        parent_res : NodeResults
            The results of the tests performed on the parent node
        node_res : NodeResults
            The results of the tests performed on the current node
        num_samps : int
            The number of samples tested
        num_tests : The number of haplotypes tested

        Returns
        -------
        bool
            True if the branch should be terminated, False otherwise
        """
        t_stat = None
        if parent_res:
            # before we do any calculations, check whether the effect sizes have
            # improved and return True if they haven't
            if (
                np.isnan(node_res.beta)
                or (np.abs(node_res.beta) - np.abs(parent_res.beta)) >= 0
            ):
                # terminate if the effect sizes have gone in the opposite direction
                self.log.debug("Terminated b/c effect size went in wrong direction")
                return True
            # perform a two tailed, two-sample t-test using the difference of the effect sizes
            # first, we compute the standard error of the difference of the effect sizes
            std_err = np.sqrt(((node_res.stderr ** 2) + (parent_res.stderr ** 2)) / 2)
            if std_err == 0:
                # if we have a standard error of 0, then we already know the result is
                # significant! It doesn't matter what the effect sizes are b/c t_stat
                # will be inf
                return False
            # then, we compute the test statistic
            # use np.abs to account for the directions that the effect size may take
            t_stat = np.abs(node_res.beta - parent_res.beta) / std_err
            # use a one-tailed test here b/c either the effect size becomes more
            # negative or it becomes more positive
            pval = t_dist.cdf(-t_stat, df=2 * (num_samps - 2))
        else:
            # parent_res = None when the parent node is the root node
            pval = node_res.pval
        if np.isnan(pval):
            raise ValueError(
                "Encountered an nan p-value! Check your data for irregularities."
            )
        # correct for multiple hypothesis testing
        # For now, we use the Bonferroni correction
        if pval >= (self.method.pval_thresh / num_tests):
            self.log.debug("Terminated with t-stat {} and p-value {}".format(t_stat, pval))
            return True
        self.log.debug("Significant with t-stat {} and p-value {}".format(t_stat, pval))
        return False
