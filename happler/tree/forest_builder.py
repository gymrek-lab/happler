from __future__ import annotations
from logging import Logger

import numpy as np
import numpy.typing as npt
import statsmodels.api as sm
from haptools.logging import getLogger
from haptools.data import Genotypes, Phenotypes

from .tree import Tree
from .haplotypes import Haplotypes
from .tree_builder import TreeBuilder


class ForestBuilder:
    """
    Creates a ForestBuilder object from a provided TreeBuilder

    Attributes
    ----------
    tree_builder: TreeBuilder
        A TreeBuilder instance that we can run repeatedly
    log: Logger
        A logging instance for recording debug statements.
    """

    def __init__(
        self,
        tree_builder: TreeBuilder,
        num_bins: int,
        max_iterations: int = None,
        log: Logger = None,
    ):
        self.tree_builder = tree_builder
        self.max_iterations = max_iterations
        self.trees = [None]*num_bins
        self.haplotypes = [None]*num_bins
        self.genotypes = self.tree_builder.gens
        self.phenotypes = self.tree_builder.phens
        self.log = log or getLogger(self.__class__.__name__)

    def run(self):
        """
        Build a bunch of trees
        """
        iterations = self.max_iterations
        tree_idx = 0
        while iterations:
            self.log.debug(f"Iteration {self.max_iterations - iterations}: tree {tree_idx}")
            self.tree_builder.phens = self.get_residuals(tree_idx)
            self.tree_builder.run()
            self.trees[tree_idx] = self.tree_builder.tree
            self.haplotypes[tree_idx] = Haplotypes.from_tree(
                None, self.trees[tree_idx], self.genotypes, self.log,
            )
            # increment tree_idx until it reaches the end
            tree_idx += 1
            if tree_idx >= len(self.trees):
                # TODO: implement a test to check whether to stop iterating
                if iterations is not None:
                    iterations -= 1
                tree_idx = 0
        return self.haplotypes

    def get_residuals(self, tree_idx: int) -> Phenotypes:
        """
        Retrieve the current phenotype values, computed after regressing out the
        haplotype effects from each tree except the one indicated by tree_idx

        Parameters
        ----------
        tree_idx: int
            The index of the tree whose haplotypes to exclude from consideration

        Returns
        -------
        Phenotypes
            The residuals after regressing out the effects of the haplotypes from all
            other trees
        """
        # first: count the number of effects amongst all haplotypes
        total_num_effects = sum(
            len(hps.data) if hps is not None else 0
            for hps_idx, hps in enumerate(self.haplotypes)
            if hps_idx != tree_idx
        )
        # base case: if there are no effects to regress out
        if total_num_effects == 0:
            return self.phenotypes
        num_samps = len(self.phenotypes.samples)
        effects = np.empty((num_samps, total_num_effects), dtype=self.phenotypes.data.dtype)
        effect_arr_idx = 0
        # iterate through every tree's Haplotypes obj and transform the haps into effects
        self.log.debug(f"Extracting {total_num_effects} effects from {total_num_effects} total effects")
        for hap_idx, hap in enumerate(self.haplotypes):
            if hap is None or hap_idx == tree_idx:
                continue
            new_idx = effect_arr_idx + len(hap.data)
            effects[:, effect_arr_idx:new_idx] = hap.transform(self.genotypes).data[:,:,:2].sum(axis=2)
            effect_arr_idx = new_idx
        resids = Phenotypes(fname=None, log=self.phenotypes.log)
        resids.samples = self.phenotypes.samples
        # get residuals with effects as covariates
        resids.names = self.phenotypes.names
        self.log.debug(f"Computing residuals for tree {tree_idx}")
        resids.data = sm.OLS(self.phenotypes.data, sm.add_constant(effects)).fit().resid[:, np.newaxis]
        return resids

    def __repr__(self):
        return self.dot()

    def dot(self) -> str:
        """
        Convert the trees to a representation in the dot language
        This is useful for quickly viewing the forest on the command line

        Returns
        -------
        str
            A string representing the trees

            Nodes are labeled by their variant ID and edges are labeled by their allele
        """
        trees_str = "strict digraph {\nforcelabels=true;\nrankdir=TB;\n"
        for idx, tree in enumerate(self.trees):
            if tree is None:
                continue
            trees_str += f"subgraph tree_{idx}"+" {\nlabel=\""+f"Tree {idx}"+"\";\n"
            trees_str += "\n".join(tree.dot().split("\n")[2:])
        trees_str += "}"
        return trees_str
