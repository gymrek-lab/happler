import numpy as np
import networkx as nx

from dataclasses import dataclass
from __future__ import annotations
from collections import defaultdict

from ..data import VariantType


# We declare this class to be a dataclass to automatically define __init__ and a few
# other methods. We use frozen=True to make it immutable.
@dataclass(frozen=True)
class Node:
    """
    A node within the tree

    Attributes
    ----------
    id : str
        The variant's unique ID
    idx : int
        The index of the variant within the genotype data
    pos : int
        The chromosomal position of the variant
    """

    idx: int
    id: int
    pos: int

    @classmethod
    def from_np(cls, np_mixed_arr_var: np.void, idx: int) -> Node:
        """
        Convert a numpy mixed array variant record into a Node

        Parameters
        ----------
        np_mixed_arr_var : np.void
            A numpy mixed array variant record with entries 'id' and 'pos'
        idx : int
            See :py:attr:`~.Node.idx`

        Returns
        -------
        Node
            The converted Node
        """
        return cls(idx=idx, id=np_mixed_arr_var["id"], pos=np_mixed_arr_var["pos"])

    @property
    def ID(self):
        return self.id

    @property
    def POS(self):
        return self.pos


# We declare this class to be a dataclass to automatically define __init__ and a few
# other methods. We use frozen=True to make it immutable.
@dataclass(frozen=True)
class Haplotype:
    """
    A haplotype within the tree

    Attributes
    ----------
    nodes : tuple[tuple[Node, int]]
        An ordered collection of pairs, where each pair is a node and its allele
    """

    # TODO: consider using a named tuple? ugh maybe we want both a Node class and a
    # Variant class?
    nodes: tuple[tuple[Node, int]]

    @classmethod
    def from_node(cls, node: Node, allele: int) -> Haplotype:
        """
        Create a new haplotype with a single node entry

        Parameters
        ----------
        node : Node
            The initializing node for this haplotype
        allele : int
            The allele associated with node

        Returns
        -------
        Haplotype
            The newly created haplotype object containing node and allele
        """
        return cls(((Node, allele),))

    def append(self, node: Node, allele: int) -> Haplotype:
        """
        Append a new node to this haplotype

        Parameters
        ----------
        node : Node
            The node to add to this haplotype
        allele : int
            The allele associated with this node

        Returns
        -------
        Haplotype
            A new haplotype object extended by the node and its allele
        """
        return Haplotype(self.nodes + ((Node, allele),))

    @property
    def node_indices(self) -> tuple[int]:
        """
        Get the indices of the nodes in this haplotype

        Returns
        -------
        tuple
            The indices of the nodes
        """
        return tuple(node[0].idx for node in self.nodes)


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

    def _add_root_node(self, root: Node):
        """
        Initialize the tree with a root node

        Parameters
        ----------
        root : Node
            See :py:attr:`_.Tree.root`
        """
        # TODO: consider reducing code duplication between this method and add_node()
        self.variant_locs[root].add(self.num_nodes)
        self.graph.add_node(self.num_nodes, idx=root.idx, label=root.id)

    def add_node(self, node: Node, parent_idx: int, allele: int) -> int:
        """
        Add node to the tree

        Parameters
        ----------
        node : Node
            The node to add to the tree
        parent_idx : int
            The index of the node under which to place the new node
        allele : int
            The allele for the edge from parent to node

        Returns
        -------
        int
            The index of the new node within the tree
        """
        new_node_idx = self.num_nodes
        self.graph.add_node(new_node_idx, idx=node.idx, label=node.id)
        self.variant_locs[node].add(new_node_idx)
        self.graph.add_edge(parent_idx, new_node_idx, label=allele)
        return new_node_idx

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
