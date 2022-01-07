from __future__ import annotations
from collections import defaultdict

import numpy as np
import networkx as nx

from .variant import Variant


class Tree:
    """
    A tree where

    1. nodes are variants (and also haplotypes)
    2. each node can have at most two branches for each of its alleles

    Attributes
    ----------
    graph: nx.DiGraph
        The underlying directed acyclic graph representing the tree
    variant_locs: dict
        The indices of each variant within the tree's list of nodes
    """

    def __init__(self, root: None):
        self.graph = nx.DiGraph()
        self.variant_locs = defaultdict(set)
        self._add_root_node(root)

    @property
    def num_nodes(self):
        return self.graph.number_of_nodes()

    def _add_root_node(self, root: Variant):
        """
        Initialize the tree with a root node

        Parameters
        ----------
        root : Variant
            See :py:attr:`_.Tree.root`
        """
        # TODO: consider reducing code duplication between this method and add_node()
        self.variant_locs[root].add(self.num_nodes)
        self.graph.add_node(self.num_nodes, idx=root.idx, label=root.id)

    def add_node(self, node: Variant, parent_idx: int, allele: int, pval: float) -> int:
        """
        Add node to the tree

        Parameters
        ----------
        node : Variant
            The node to add to the tree
        parent_idx : int
            The index of the node under which to place the new node
        allele : int
            The allele for the edge from parent to node
        pval : float
            The p-value of the haplotype once this variant is included

        Returns
        -------
        int
            The index of the new node within the tree
        """
        new_node_idx = self.num_nodes
        self.graph.add_node(new_node_idx, idx=node.idx, label=node.id, pval=pval)
        self.variant_locs[node].add(new_node_idx)
        self.graph.add_edge(parent_idx, new_node_idx, label=allele)
        return new_node_idx

    def haplotypes(self) -> list[list[tuple[str, float]]]:
        """
        Return the haplotypes at the leaves of this tree

        Returns
        -------
        list[list[tuple[str, float]]]
            A list of haplotypes, where each haplotype consists of a list of tuples: a variant ID and a score
        """
        pass

    def ascii(self) -> str:
        """
        Convert the tree to its ascii representation
        This is useful for quickly viewing the tree on the command line

        Returns
        -------
        str
            A string representing the tree

            Nodes are labeled by their variant ID and edges are labeled by their allele
        """
        pass
