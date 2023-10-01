from __future__ import annotations
from logging import Logger, DEBUG
from collections import defaultdict, deque

import networkx as nx
from haptools.logging import getLogger

from .variant import Variant
from .assoc_test import NodeResults


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

    def __init__(self, log: Logger = None):
        self.graph = nx.DiGraph()
        self.variant_locs = defaultdict(set)
        self.log = log or getLogger(self.__class__.__name__)
        self._add_root_node()

    def __repr__(self):
        return self.dot()

    @property
    def num_nodes(self):
        return self.graph.number_of_nodes()

    @property
    def num_variants(self):
        return len(self.variant_locs)

    def _add_root_node(self):
        """
        Initialize the tree with a special root node
        Note that this node does not represent a variant of any kind
        """
        self.graph.add_node(0, variant=None, label=None, allele=None, results=None)

    def add_node(
        self, node: Variant, parent_idx: int, allele: int, results: NodeResults = None
    ) -> int:
        """
        Add node to the tree

        Parameters
        ----------
        node : Variant
            The node to add to the tree
        parent_idx : int
            The index of the node under which to place the new node
        allele : int
            The allele for the edge from a parent node to this node
        results : np.void, optional
            The results (beta, pval) of the association test for this haplotype once
            this variant is included

        Returns
        -------
        int
            The index of the new node within the tree
        """
        if self.graph.out_degree(parent_idx) >= 2:
            raise ValueError(
                "This parent node has reached its capacity already. It cannot take"
                " more children."
            )
        new_node_idx = self.num_nodes
        label = ""
        if node is not None:
            label = node.id
        self.graph.add_node(
            new_node_idx, variant=node, label=label, allele=allele, results=results
        )
        self.variant_locs[node].add(new_node_idx)
        self.graph.add_edge(parent_idx, new_node_idx, label=allele)
        return new_node_idx

    def remove_leaf_node(self, node_idx: int):
        """
        Remove a leaf node from the tree

        Parameters
        ----------
        node_idx : int
            The index of the node to remove
        """
        # check that this node is a leaf
        if self.graph.out_degree[node_idx]:
            raise ValueError("Cannot remove non-leaf node.")
        variant = self.graph.nodes[node_idx]["variant"]
        # remove the node from the graph
        self.graph.remove_node(node_idx)
        self.variant_locs[variant].remove(node_idx)

    def siblings(self, node_idx: int) -> list[Variant]:
        """
        Locate sibling(s) of this node in the tree

        Parameters
        ----------
        node_idx : int
            The index of a node in the tree

        Returns
        -------
        list[Variant]
            The variants at the sibling nodes
        """
        try:
            parent = next(self.graph.predecessors(node_idx))
        except StopIteration:
            return {}
        return {
            node: self.graph.nodes[node]
            for node in self.graph.successors(parent)
            if node != node_idx
        }

    def leaves(self, from_root: bool = False):
        """
        Return all leaves of this tree

        Parameters
        ----------
        from_root: bool, optional
            Whether to only return leaves attached to the root of the tree

        Returns
        -------
        tuple[int, Variant, int]
            The variant at each leaf node of the tree. Returns the index in the tree,
            the Variant object, and the allele of the variant.
        """
        if from_root:
            from_root = lambda node: next(self.graph.predecessors(node)) == 0
        else:
            from_root = lambda node: True
        return {
            node: self.graph.nodes[node]
            for node, degree in self.graph.out_degree
            if degree == 0 and node != 0 and from_root(node)
        }

    def haplotypes(self, root: int = 0) -> list[deque[dict]]:
        """
        Return the haplotypes at the leaves of the tree rooted at the index "root"

        Returns
        -------
        list[deque[dict]]
            A list of haplotypes, where each haplotype consists of a list of dictionaries.

            Each dictionary contains the node and all of its attributes.
        """
        # how many children does this node have?
        num_children = self.graph.out_degree(root)
        # check: is the root index 0 or some other int?
        if root:
            root_node = deque([self.graph.nodes[root]])
        elif num_children:
            # this root is actually the absolute root, which doesn't represent a real
            # variant, so we just use an empty deque in that case
            root_node = deque([])
        else:
            # if the root is actually the absolute root and there aren't any children,
            # just return an empty list
            return []
        # first, check that this node is not a leaf
        if num_children:
            return [
                root_node + path
                for child in self.graph.successors(root)
                for path in self.haplotypes(child)
            ]
        # if it's a leaf, just output a single node
        return [root_node]

    def dot(self) -> str:
        """
        Convert the tree to its representation in the dot language
        This is useful for quickly viewing the tree on the command line

        Returns
        -------
        str
            A string representing the tree

            Nodes are labeled by their variant ID and edges are labeled by their allele
        """
        dot = nx.drawing.nx_pydot.to_pydot(self.graph)
        dot.obj_dict["attributes"]["forcelabels"] = "true"
        # iterate through all of the nodes, treating the root specially
        for idx, node in enumerate(dot.get_nodes()):
            # node.set_name(node.get('label'))
            attrs = node.get_attributes()
            # check: does this node have a valid variant attached to it?
            if attrs["variant"] == "None" and idx == 0:
                # treat the root node specially, since it isn't a real variant
                node.obj_dict["attributes"] = {"label": "root"}
            elif self.log.getEffectiveLevel() == DEBUG:
                node.obj_dict["attributes"] = {
                    "label": attrs["label"] + "\n" + attrs["results"]
                }
            else:
                node.obj_dict["attributes"] = {"label": attrs["label"]}

        return dot.to_string()
