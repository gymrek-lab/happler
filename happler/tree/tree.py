import numpy as np
import networkx as nx

from ..data import VariantType


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
    type : VariantType
        The type of variant
    _tree_idx : int
        The index of the variant within the tree
    """

    __slots__ = "idx", "id", "pos", "type"

    def __init__(
        self,
        idx: int,
        ID: str,
        pos: int,
        variant_type: VariantType = VariantType("snp"),
    ):
        self.idx = idx
        self.id = ID
        self.pos = pos
        self.type = variant_type
        self._tree_idx = None

    @classmethod
    def from_np(cls, idx: int, np_mixed_arr_var: np.void):
        """
        Convert a numpy mixed array variant record into a Node

        Parameters
        ----------
        idx : int
            See :py:attr:`~.Node.idx`
        np_mixed_arr_var : np.void
            A numpy mixed array variant record with entries 'id' and 'pos'

        Returns
        -------
        Node
            The converted Node
        """
        return cls(idx=idx, ID=np_mixed_arr_var["id"], pos=np_mixed_arr_var["pos"])

    @property
    def ID(self):
        return self.id

    @property
    def POS(self):
        return self.pos

    def __repr__(self):
        return self.id

    def __hash__(self):
        return self._tree_idx

    def __eq__(self, other):
        return self.idx == other.idx


class Tree:
    """
    A tree where

    1. nodes are variants (and also haplotypes)
    2. each node can have at most two branches for each of its alleles

    Attributes
    ----------
    graph: nx.DiGraph
        The underlying directed acyclic graph representing the tree
    """

    def __init__(self):
        self.graph = nx.DiGraph()
