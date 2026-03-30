from __future__ import annotations
import math
import logging
from logging import Logger

import numpy as np
import numpy.typing as npt
from haptools.logging import getLogger
from haptools.data import Genotypes, Phenotypes
from haptools.ld import pearson_corr_ld

from .tree import Tree
from .variant import Variant
from .haplotypes import Haplotype
from .terminator import Terminator, BICTerminator, TTestTerminator
from .assoc_test import (
    AssocTest,
    NodeResults,
    NodeResultsBIC,
    NodeResultsExtra,
    AssocTestSimpleFastBIC,
    AssocTestSimpleSMTScore,
    AssocTestSimpleCovariates,
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
        method: AssocTest = AssocTestSimpleFastBIC(),
        terminator: Terminator = BICTerminator(),
        indep_thresh: float = 15,
        ld_prune_thresh: float = None,
        covariance_correction: float = True,
        log: Logger = None,
    ):
        self.gens = genotypes
        self.phens = phenotypes
        self.maf = maf
        self.method = method
        self.terminator = terminator
        self.indep_thresh = indep_thresh
        self.results_type = method.results_type
        self.tree = None
        self._split_method = self._find_split_rigid
        split_method = "rigid"
        self.ld_prune_thresh = ld_prune_thresh
        self.ranking_val = (
            "bic" if isinstance(self.method, AssocTestSimpleFastBIC) else "pval"
        )
        # for now, let's comment this out because we want to try the rigid strategy
        # if self.ld_prune_thresh is not None:
        #     self._split_method = self._find_split_flexible
        #     split_method = "flexible"
        self.covariance_correction = covariance_correction
        self.log = log or getLogger(self.__class__.__name__)
        self.log.info(f"Using {split_method} branching strategy")

    def __repr__(self):
        return str((self.gens, self.phens))

    @staticmethod
    def _make_align_mask(
        num_variants: int, valid_idxs: npt.NDArray[np.int_]
    ) -> npt.NDArray[np.bool_]:
        """
        Build a boolean mask of length num_variants indicating which global variant indices
        are represented by the columns of the current X matrix passed to AssocTest.run().
        """
        align_mask = np.zeros(num_variants, dtype=bool)
        align_mask[valid_idxs] = True
        return align_mask

    def run(self):
        """
        Run the tree builder and create a tree rooted at the provided variant
        """
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
        parent_maf_mask: dict = None,
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
        parent_maf_mask : dict, optional
            A dict mapping allele (0 or 1) to an array of valid variant indices
            (into the full genotype matrix) that passed MAF filtering at the parent
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
        vals, child_maf_mask = self._split_method(parent_hap, parent_res, parent_maf_mask)
        if vals is not None:
            for variant, allele, results in vals:
                if variant is None:
                    # there were no significant variants!
                    continue
                new_node_idx = self.tree.add_node(variant, parent_idx, allele, results)
                if self.log.getEffectiveLevel() == logging.DEBUG:
                    self.log.debug(self.tree.dot())
                # create a new Haplotype with the variant-allele pair added
                variant_gts = self.gens.data[:, variant.idx, :2] == allele
                new_parent_hap = parent_hap.append(variant, allele, variant_gts)
                self._create_tree(new_parent_hap, new_node_idx, results, child_maf_mask)

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
                sib_p = sibling["results"].bic
                leaf_p = leaf["results"].bic
                if sib_p > leaf_p and not math.isclose(sib_p, leaf_p):
                    self.log.debug(
                        f"Left leaf {leaf_var.id} unpruned since it has a better bic"
                    )
                    # keep it if our value is better
                    continue
                elif math.isclose(sib_p, leaf_p) and leaf["allele"] == 1:
                    self.log.debug(
                        f"Left leaf {leaf_var.id} unpruned since it's for the ALT allele"
                    )
                    # also if the values are the same but our effect size is positive
                    continue
            # step 3: get the genotypes for the leaf node and its sibling
            leaf_gts = self.gens.data[:, leaf_var.idx, :] == leaf["allele"]
            sibling_gts = (
                self.gens.data[:, sibling["variant"].idx, :] == sibling["allele"]
            )
            # step 4: check whether the leaf node is in strong LD with its sibling
            ld = pearson_corr_ld(leaf_gts.sum(axis=1), sibling_gts.sum(axis=1)) ** 2
            if ld > self.ld_prune_thresh:
                self.log.debug(f"Pruning {leaf_var.id} with LD {ld}")
                count += 1
                # step 5: if it is, remove it
                self.tree.remove_leaf_node(leaf_idx)
            else:
                self.log.debug(f"Left leaf {leaf_var.id} (with LD {ld}) unpruned")
        self.log.debug(
            f"Pruned {count} leaves with LD > {self.ld_prune_thresh} with their siblings"
        )

    def _find_split_flexible(
        self,
        parent: Haplotype,
        parent_res: NodeResults = None,
        parent_maf_mask: dict = None,
    ) -> tuple[Variant, np.void]:
        """
        Find the variant/allele that best fits under the parent_idx node

        Parameters
        ----------
        parent : Haplotype
            The haplotype containing all variants up to (and including) the parent
        parent_res : NodeResults, optional
            The results of the tests performed on the parent node
        parent_maf_mask : dict, optional
            A dict mapping allele (0 or 1) to an array of valid variant indices
            (into the full genotype matrix) that passed MAF filtering at the parent

        Returns
        -------
        tuple[Variant, float]
            The variant that best fits under the parent node with the allele edge AND
            the results (ex: beta, bic) of the haplotype association test after
            incorporating that variant
        """
        final_to_return = []
        num_samps = len(self.gens.samples)
        num_variants = self.gens.data.shape[1]
        parent_indices_set = set(parent.node_indices)
        child_maf_mask = {}
        # iterate through the two possible alleles and try all SNPs with that allele
        alleles = (0, 1)
        for allele in alleles:
            # step 0: determine variant indices to consider for this allele; start from
            # parent mask and exclude parent nodes
            if parent_maf_mask is not None:
                valid = parent_maf_mask[allele]
                if parent_indices_set:
                    valid = valid[~np.isin(valid, list(parent_indices_set))]
            else:
                all_indices = np.arange(num_variants)
                if parent_indices_set:
                    valid = all_indices[~np.isin(all_indices, list(parent_indices_set))]
                else:
                    valid = all_indices
            # step 1: transform the GT matrix and sum along ploidy axis using JAX JIT
            hap_mat_sum = parent.transform_and_sum(
                self.gens, allele, idxs=valid, remove_self=False
            )
            # step 2: exclude any haplotypes that are now too rare
            if self.maf is not None:
                ref_af = hap_mat_sum.sum(axis=0) / hap_mat_sum.shape[0] / 2
                maf = np.minimum(ref_af, 1 - ref_af)
                maf_pass = np.nonzero(maf >= self.maf)[0]
                child_maf_mask[allele] = valid[maf_pass]
                if len(maf_pass) < hap_mat_sum.shape[1]:
                    self.log.debug(
                        f"Considering {len(maf_pass)} variants for allele {allele}"
                    )
                    hap_mat_sum = hap_mat_sum[:, maf_pass]
            else:
                child_maf_mask[allele] = valid
            if hap_mat_sum.shape[1] == 0:
                # if there weren't any genotypes left, just return None
                self.log.debug(f"No variants passed --hap-maf for allele {allele}")
                final_to_return.append((None, allele, None))
                continue
            parent_corr = None
            align_mask = self._make_align_mask(num_variants, child_maf_mask[allele])
            if isinstance(self.method, AssocTestSimpleSMTScore) and not (
                parent_res is None
            ):
                if self.covariance_correction:
                    parent_corr = pearson_corr_ld(hap_mat_sum, parent.data.sum(axis=1))
                # step 3: run all association tests on all of the haplotypes
                results = self.method.run(
                    hap_mat_sum,
                    self.phens.data[:, 0],
                    parent_res=parent_res,
                    parent_corr=parent_corr,
                    align_to_mask=align_mask,
                )
                # step 4: record the best t-score among all the SNPs with this allele
                best_var_idx = int(results.data["tscore"].argmax())
            else:
                # step 3: run all association tests on all of the haplotypes
                results = self.method.run(
                    hap_mat_sum,
                    self.phens.data[:, 0],
                    align_to_mask=align_mask,
                )
                # step 4: record the best p-value among all the SNPs with this allele
                best_var_idx = int(results.data[self.ranking_val].argmin())
            node_res = self.results_type.from_np(results.data[best_var_idx])
            num_tests = len(parent.nodes) + 1
            # step 5: retrieve the Variant with the best value
            best_variant = Variant.from_np(self.gens.variants[best_var_idx], best_var_idx)
            self.log.debug("Chose variant {}".format(best_variant.id))
            # step 6: check whether we don't get a stronger effect by treating this variant
            # as independently causal
            if parent_res is not None:
                allele_gts = (self.gens.data[:, best_var_idx] == allele).sum(axis=1)[
                    :, np.newaxis
                ]
                # y ~ h_parent + z_child VS y ~ h_hap
                hap_indep_effect = NodeResultsExtra.from_np(
                    AssocTestSimpleCovariates(covars=allele_gts, with_bic=True)
                    .run(
                        parent.data.sum(axis=1)[:, np.newaxis],
                        self.phens.data[:, 0],
                    )
                    .data[0]
                )
                if BICTerminator(bf_thresh=self.indep_thresh).check(
                    hap_indep_effect,
                    node_res,
                    results,
                    best_var_idx,
                    num_samps,
                    num_tests,
                ):
                    self.log.debug(
                        "Terminating because the hap had a BIC too similar to one with "
                        "just the parent + child"
                    )
                    final_to_return.append((None, allele, node_res))
                    continue
                self.log.debug(
                    "The haplotype had a much better BIC than in an additive model "
                    "with the allele and parent"
                )
            # step 7: check if this allele is significant and whether we should terminate the branch
            self.log.debug(
                "Testing variant {} / allele {} with parent_res {} and node_res {}".format(
                    best_variant.id, allele, parent_res, node_res
                )
            )
            if self.terminator.check(
                parent_res,
                node_res,
                results,
                best_var_idx,
                num_samps,
                num_tests,
            ):
                final_to_return.append((None, allele, node_res))
                continue
            final_to_return.append((best_variant, allele, node_res))
        return final_to_return, child_maf_mask

    def _find_split_rigid(
        self,
        parent: Haplotype,
        parent_res: NodeResults = None,
        parent_maf_mask: dict = None,
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
        parent_maf_mask : dict, optional
            A dict mapping allele (0 or 1) to an array of valid variant indices
            (into the full genotype matrix) that passed MAF filtering at the parent

        Returns
        -------
        tuple[Variant, float]
            The variant that best fits under the parent node with the allele edge AND
            the results (ex: beta, bic) of the haplotype association test after
            incorporating that variant
        """
        final_to_return = []
        num_samps = len(self.gens.samples)
        num_variants = self.gens.data.shape[1]
        parent_indices_set = set(parent.node_indices)
        results = {}
        best_var_by_allele = {}
        maf_mask = {}
        parent_corr = {}
        # iterate through the two possible alleles and try all SNPs with that allele
        alleles = (0, 1)
        for allele in alleles:
            # step 0: determine variant indices to consider for this allele; start from
            # parent mask and exclude parent nodes
            if parent_maf_mask is not None:
                valid = parent_maf_mask[allele]
                if parent_indices_set:
                    valid = valid[~np.isin(valid, list(parent_indices_set))]
            else:
                all_indices = np.arange(num_variants)
                if parent_indices_set:
                    valid = all_indices[~np.isin(all_indices, list(parent_indices_set))]
                else:
                    valid = all_indices
            # step 1: transform the GT matrix and sum along ploidy axis using JAX JIT
            hap_mat_sum = parent.transform_and_sum(
                self.gens, allele, idxs=valid, remove_self=False
            )
            # step 2: exclude any haplotypes that are too rare
            if self.maf is not None:
                ref_af = hap_mat_sum.sum(axis=0) / hap_mat_sum.shape[0] / 2
                maf = np.minimum(ref_af, 1 - ref_af)
                maf_pass = np.nonzero(maf >= self.maf)[0]
                maf_mask[allele] = valid[maf_pass]
                if len(maf_pass) < hap_mat_sum.shape[1]:
                    self.log.debug(
                        f"Considering {len(maf_pass)} variants for allele {allele}"
                    )
                    hap_mat_sum = hap_mat_sum[:, maf_pass]
            else:
                maf_mask[allele] = valid
            if hap_mat_sum.shape[1] == 0:
                # if there weren't any genotypes left, just return None
                self.log.debug(f"No variants passed --hap-maf for allele {allele}")
                final_to_return.append((None, allele, None))
                continue
            parent_corr[allele] = None
            align_mask = self._make_align_mask(num_variants, maf_mask[allele])
            if isinstance(self.method, AssocTestSimpleSMTScore) and not (
                parent_res is None
            ):
                if self.covariance_correction:
                    parent_corr[allele] = pearson_corr_ld(
                        hap_mat_sum, parent.data.sum(axis=1)
                    )
                # step 3: run all association tests on all of the haplotypes
                results[allele] = self.method.run(
                    hap_mat_sum,
                    self.phens.data[:, 0],
                    parent_res=parent_res,
                    parent_corr=parent_corr[allele],
                    align_to_mask=align_mask,
                )
                # step 4: record the best t-score among all the SNPs with this allele
                best_var_by_allele[allele] = int(results[allele].data["tscore"].argmax())
            else:
                # step 3: run all association tests on all of the haplotypes
                results[allele] = self.method.run(
                    hap_mat_sum,
                    self.phens.data[:, 0],
                    align_to_mask=align_mask,
                )
                # step 4: record the best t-score among all the SNPs with this allele
                best_var_by_allele[allele] = int(
                    results[allele].data[self.ranking_val].argmin()
                )
        # exit if neither of the alleles worked
        if not len(best_var_by_allele):
            return None, maf_mask
        # step 5: find the index of the best variant within the haplotype matrix
        if isinstance(self.method, AssocTestSimpleSMTScore):
            best_allele = max(
                best_var_by_allele,
                key=lambda a: results[a].data["tscore"][best_var_by_allele[a]],
            )
        else:
            best_allele = min(
                best_var_by_allele,
                key=lambda a: results[a].data[self.ranking_val][best_var_by_allele[a]],
            )
        best_var_idx = best_var_by_allele[best_allele]
        # For rigid mode: find the same variant in the other allele's results
        other_allele = int(not best_allele)
        # With global-aligned results, both alleles can be handled using the same global best_var_idx.
        # If a SNP wasn't tested for an allele (e.g. no variants passed MAF), skip that allele.
        best_res_idx = {best_allele: best_var_idx}
        if other_allele in results:
            best_res_idx[other_allele] = best_var_idx
        else:
            final_to_return.append((None, allele, None))
        num_tests = len(parent.nodes) + 1
        # step 6: retrieve the Variant with the best value
        best_variant = Variant.from_np(self.gens.variants[best_var_idx], best_var_idx)
        self.log.debug("Chose variant {}".format(best_variant.id))
        # step 6: check the MAFs of the haplotypes we created
        # if the best variant was filtered out for this allele due to low MAF
        # iterate through all of the alleles of the best variant and check if they're
        # significant
        for allele in best_res_idx:
            best_results = results[allele].data[best_var_idx]
            node_res = self.results_type.from_np(best_results)
            # step 7: check whether we don't get a stronger effect by treating this variant
            # as independently causal
            if parent_res is not None:
                allele_gts = self.gens.data[:, best_var_idx].sum(axis=1)[:, np.newaxis]
                # y ~ h_parent + z_child VS y ~ h_hap
                # delta BIC = first minus second
                hap_indep_effect = NodeResultsExtra.from_np(
                    AssocTestSimpleCovariates(covars=allele_gts, with_bic=True)
                    .run(
                        parent.data.sum(axis=1)[:, np.newaxis],
                        self.phens.data[:, 0],
                    )
                    .data[0]
                )
                if BICTerminator(bf_thresh=self.indep_thresh).check(
                    hap_indep_effect,
                    node_res,
                    results[allele],
                    best_var_idx,
                    num_samps,
                    num_tests,
                ):
                    self.log.debug(
                        "Terminating because the hap had a BIC too similar to one with "
                        f"just the parent + child for allele {allele}"
                    )
                    final_to_return.append((None, allele, node_res))
                    continue
                self.log.debug(
                    "The haplotype had a much better BIC than in an additive model "
                    f"with the allele and parent for allele {allele}"
                )
            # step 8: check whether we should terminate the branch
            self.log.debug(
                "Testing variant {} / allele {} with parent_res {} and node_res {}".format(
                    best_variant.id, allele, parent_res, node_res
                )
            )
            if self.terminator.check(
                parent_res,
                node_res,
                results[allele],
                best_var_idx,
                num_samps,
                num_tests,
            ):
                final_to_return.append((None, allele, node_res))
                continue
            final_to_return.append((best_variant, allele, node_res))
        return final_to_return, maf_mask
