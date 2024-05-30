from __future__ import annotations
from logging import Logger

import numpy as np
from haptools.logging import getLogger
from haptools.data import Genotypes, Phenotypes
from haptools.ld import pearson_corr_ld

from .tree import Tree
from .variant import Variant
from .haplotypes import Haplotype
from .terminator import Terminator, TTestTerminator
from .assoc_test import (
    AssocTest,
    NodeResults,
    NodeResultsExtra,
    NodeResultsExtraTScore,
    NodeResultsTScore,
    AssocTestSimple,
    AssocTestSimpleSMTScore,
)


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
    ld_prune_thresh: float, optional
        Any leaf nodes with a greater LD with their sibling than this value will be pruned
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
        maf: float = None,
        method: AssocTest = AssocTestSimple(),
        terminator: Terminator = TTestTerminator(),
        ld_prune_thresh: float = None,
        log: Logger = None,
    ):
        self.gens = genotypes
        self.phens = phenotypes
        self.maf = maf
        self.method = method
        self.terminator = terminator
        self.results_type = method.results_type
        self.tree = None
        self._split_method = self._find_split_rigid
        split_method = "rigid"
        self.ld_prune_thresh = ld_prune_thresh
        # for now, let's comment this out because we want to try the rigid strat
        # if self.ld_prune_thresh is not None:
            # self._split_method = self._find_split_flexible
            # split_method = "flexible"
        self.log = log or getLogger(self.__class__.__name__)
        self.log.info(f"Using {split_method} branching strategy")

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
        self.tree = Tree(log=self.log)
        # step two: create a haplotype
        # we're at the root of tree, so we need to create a new, empty haplotype
        parent_hap = Haplotype(num_samples=len(self.gens.samples))
        # step three: create the rest of the tree
        self._create_tree(parent_hap, parent_idx=0)
        # step four: prune nodes from the tree that are in strong LD
        self.prune_tree()
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
        for variant, allele, results in self._split_method(parent_hap, parent_res):
            if variant is None:
                # there were no significant variants!
                continue
            new_node_idx = self.tree.add_node(variant, parent_idx, allele, results)
            self.log.debug(self.tree.dot())
            # create a new Haplotype with the variant-allele pair added
            variant_gts = self.gens.data[:, variant.idx, :2] == allele
            new_parent_hap = parent_hap.append(variant, allele, variant_gts)
            self._create_tree(new_parent_hap, new_node_idx, results)

    def prune_tree(self, from_root: bool = True):
        """
        Remove any leaf nodes that are in strong LD with their sibling branches

        Parameters
        ----------
        from_root: bool, optional
            Whether to only prune leaves attached to the root of the tree
        """
        if self.ld_prune_thresh is None:
            return
        count = 0
        leaves = self.tree.leaves(from_root=from_root)
        self.log.debug(f"Considering {len(leaves)} leaves for pruning")
        # step 1: get all leaf nodes and their siblings
        for leaf_idx, leaf in leaves.items():
            leaf_var = leaf["variant"]
            # step 2: get the siblings of this leaf and check if they are leaves
            sibs = list(self.tree.siblings(leaf_idx).items())
            if not len(sibs):
                continue
            # just use the first sibling for now
            # TODO: use a for-loop if we allow more than two branches per node
            sib_idx, sibling = sibs[0]
            if sib_idx in leaves:
                sib_p = sibling["results"].pval
                leaf_p = leaf["results"].pval
                if sib_p > leaf_p and not np.isclose(sib_p, leaf_p):
                    self.log.debug(
                        f"Left leaf {leaf_var.id} unpruned since it has a better pval"
                    )
                    # keep it if our p-value is better
                    continue
                elif np.isclose(sib_p, leaf_p) and leaf["results"].beta > 0:
                    self.log.debug(
                        f"Left leaf {leaf_var.id} unpruned since it's beta is positive"
                    )
                    # also if the p-values are the same but our effect size is positive
                    continue
            # step 3: get the genotypes for the leaf node and its sibling
            leaf_gts = self.gens.data[:, leaf_var.idx, :] == leaf["allele"]
            sibling_gts = (
                self.gens.data[:, sibling["variant"].idx, :] == sibling["allele"]
            )
            # step 4: check whether the leaf node is in strong LD with its sibling
            ld = pearson_corr_ld(leaf_gts.sum(axis=1), sibling_gts.sum(axis=1))
            if np.abs(ld) > self.ld_prune_thresh:
                self.log.debug(f"Pruning {leaf_var.id} with LD {ld}")
                count += 1
                # step 5: if it is, remove it
                self.tree.remove_leaf_node(leaf_idx)
            else:
                self.log.debug(f"Left leaf {leaf_var.id} (with LD {ld}) unpruned")
        self.log.info(
            f"Pruned {count} leaves with LD > {self.ld_prune_thresh} with their"
            " siblings"
        )

    def maf_mask(
        self,
        hap_matrix: npt.NDArray[np.bool_],
    ) -> npt.NDArray:
        """
        Check the minor allele frequency of each haplotype in the provided hap matrix

        Parameters
        ----------
        hap_matrix: npt.NDArray[np.bool_]
            An array of haplotype "genotypes" of shape: num_samples x num_haplotypes

        Returns
        -------
            An integer mask denoting the indices of the haplotypes that passed the
            MAF threshold
        """
        if self.maf is None:
            return np.arange(hap_matrix.shape[1])
        num_strands = 2 * hap_matrix.shape[0]
        # TODO: make this work for multi-allelic variants, too?
        ref_af = hap_matrix.sum(axis=(0, 2)) / num_strands
        maf = np.array([ref_af, 1 - ref_af]).min(axis=0)
        common_variants = maf >= self.maf
        return np.nonzero(common_variants)[0]

    def _find_split_flexible(
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
        # iterate through the two possible alleles and try all SNPs with that allele
        alleles = (0, 1)
        for allele in alleles:
            # step 1: transform the GT matrix into a matrix of common haplotypes
            hap_matrix = parent.transform(self.gens, allele)
            maf_mask = self.maf_mask(hap_matrix)
            if len(maf_mask) != hap_matrix.shape[1]:
                self.log.debug(
                    f"Considering {len(maf_mask)} variants for allele {allele}"
                )
            if hap_matrix[:, maf_mask].shape[1] == 0:
                # if there weren't any genotypes left, just return None
                yield None, allele, None
                continue
            # step 2: run all association tests on all of the haplotypes
            if isinstance(self.method, AssocTestSimpleSMTScore) and not (
                parent_res is None
            ):
                results = self.method.run(
                    hap_matrix[:, maf_mask].sum(axis=2),
                    self.phens.data[:, 0],
                    parent_res=parent_res,
                )
                # step 3: record the best t-score among all the SNPs with this allele
                best_var_idx = results.data["tscore"].argmax()
            else:
                results = self.method.run(
                    hap_matrix[:, maf_mask].sum(axis=2),
                    self.phens.data[:, 0],
                )
                # step 3: record the best p-value among all the SNPs with this allele
                best_var_idx = results.data["pval"].argmin()
            node_res = self.results_type.from_np(results.data[best_var_idx])
            best_res_idx = best_var_idx
            num_tests = len(parent.nodes) + 1
            # step 4: find the index of the best variant within the genotype matrix
            # We need to account for the rare variants that were masked out and indices
            # that we removed when running transform()
            best_var_idx += (maf_mask[best_res_idx] - len(maf_mask[:best_res_idx]))
            # There might be a faster way of doing this but for now we're just going to
            # live with it
            for gt_idx in sorted(parent.node_indices):
                # add each idx back in, so long as they are less than the target idx
                if gt_idx > best_var_idx:
                    break
                best_var_idx += 1
            # step 5: retrieve the Variant with the best p-value
            best_variant = Variant.from_np(
                self.gens.variants[best_var_idx], best_var_idx
            )
            self.log.debug("Chose variant {}".format(best_variant.id))
            # step 6: check if this allele is significant and whether we should terminate the branch
            self.log.debug(
                "Testing variant {} / allele {} with parent_res {} and node_res {}"
                .format(best_variant.id, allele, parent_res, node_res)
            )
            if self.terminator.check(parent_res, node_res, results, best_res_idx, num_samps, num_tests):
                yield None, allele, node_res
                continue
            yield best_variant, allele, node_res

    def _find_split_rigid(
        self, parent: Haplotype, parent_res: NodeResults = None
    ) -> tuple[Variant, np.void]:
        """
        Find the variant/allele that best fits under the parent_idx node

        Unlike, :py:meth:`~.TreeBuilder._find_split_flexible()`, this method
        artificially restricts the two branches from a node to representing an allele
        from the SAME SNP, instead of potentially different SNPs.

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
        maf_mask = {}
        # iterate through the two possible alleles and try all SNPs with that allele
        alleles = (0, 1)
        for allele in alleles:
            # step 1: transform the GT matrix into a matrix of common haplotypes
            hap_matrix = parent.transform(self.gens, allele)
            maf_mask[allele] = self.maf_mask(hap_matrix)
            if len(maf_mask[allele]) != hap_matrix.shape[1]:
                self.log.debug(
                    f"Considering {len(maf_mask[allele])} variants for allele {allele}"
                )
            if hap_matrix[:, maf_mask[allele]].shape[1] == 0:
                # if there weren't any genotypes left, just return None
                yield None, allele, None
                continue
            # step 2: run all association tests on all of the haplotypes
            if isinstance(self.method, AssocTestSimpleSMTScore) and not (
                parent_res is None
            ):
                results[allele] = self.method.run(
                    hap_matrix[:, maf_mask[allele]].sum(axis=2),
                    self.phens.data[:, 0],
                    parent_res=parent_res,
                )
                # also, record the best t-score among all the SNPs with this allele
                best_p_idx[allele] = results[allele].data["tscore"].argmax()
            else:
                results[allele] = self.method.run(
                    hap_matrix[:, maf_mask[allele]].sum(axis=2),
                    self.phens.data[:, 0],
                )
                # also, record the best p-value among all the SNPs with this allele
                best_p_idx[allele] = results[allele].data["pval"].argmin()
        # exit if neither of the alleles worked
        if not len(best_p_idx):
            return
        # step 3: find the index of the best variant within the haplotype matrix
        if isinstance(self.method, AssocTestSimpleSMTScore):
            best_allele = max(
                best_p_idx, key=lambda a: results[a].data["tscore"][best_p_idx[a]]
            )
        else:
            best_allele = min(
                best_p_idx, key=lambda a: results[a].data["pval"][best_p_idx[a]]
            )
        best_var_idx = best_p_idx[best_allele]
        best_res_idx = best_var_idx
        num_tests = len(parent.nodes) + 1
        # step 4: find the index of the best variant within the genotype matrix
        # We need to account for the rare variants that were masked out and indices
        # that we removed when running transform()
        best_var_idx += maf_mask[best_allele][best_res_idx] - len(maf_mask[best_allele][:best_res_idx])
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
        for allele in results:
            if allele == best_allele:
                best_allele_idx = best_p_idx[best_allele]
            else:
                best_allele_idx = np.searchsorted(maf_mask[allele], maf_mask[best_allele][best_p_idx[best_allele]])
                # if the best variant was filtered out for this allele due to low MAF
                if best_allele_idx >= len(results[allele].data):
                    yield None, allele, None
                    continue
            best_results = results[allele].data[best_allele_idx]
            node_res = self.results_type.from_np(best_results)
            # step 6: check whether we should terminate the branch
            self.log.debug(
                "Testing variant {} / allele {} with parent_res {} and node_res {}"
                .format(best_variant.id, allele, parent_res, node_res)
            )
            if self.terminator.check(parent_res, node_res, results[allele], best_res_idx, num_samps, num_tests):
                yield None, allele, node_res
                continue
            yield best_variant, allele, node_res
