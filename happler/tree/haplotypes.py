from __future__ import annotations
import numpy as np
from pathlib import Path

from ..data import Genotypes


class VariantType:
    """
    A class denoting the type of variant

    Attributes
    ----------
    type : str, optional
        The type of variant (ex: SNP, STR, etc)

        Defaults to a single nucleotide polymorphism

    """

    def __init__(self, variant_type: str = "snp"):
        # TODO: add support for STRs, SVs, etc
        supported_types = ["snp"]
        # a VariantType has a dtype of str
        self.dtype = np.dtype("S5")
        # store the type if it is supported
        variant_type = variant_type.lower()
        if variant_type in supported_types:
            self.type = variant_type
        else:
            raise ValueError("{}s are not yet supported.".format(variant_type.upper()))

    def __repr__(self):
        return self.type.upper()

    def __eq__(self, other):
        return self.type == other.type


# We declare this class to be a dataclass to automatically define __init__ and a few
# other methods. We use frozen=True to make it immutable.
@dataclass(frozen=True)
class Variant:
    """
    A variant within the genotypes matrix

    Attributes
    ----------
    id : str
        The variant's unique ID
    idx : int
        The index of the variant within the genotype data
    pos : int
        The chromosomal start position of the variant
    """

    idx: int
    id: int
    pos: int

    @classmethod
    def from_np(cls, np_mixed_arr_var: np.void, idx: int) -> Variant:
        """
        Convert a numpy mixed array variant record into a Variant

        Parameters
        ----------
        np_mixed_arr_var : np.void
            A numpy mixed array variant record with entries 'id' and 'pos'
        idx : int
            See :py:attr:`~.Variant.idx`

        Returns
        -------
        Variant
            The converted Variant
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
    nodes : tuple[tuple[Variant, int]]
        An ordered collection of pairs, where each pair is a node and its allele
    """

    # TODO: consider using a named tuple?
    nodes: tuple[tuple[Variant, int]]

    @classmethod
    def from_node(cls, node: Variant, allele: int) -> Haplotype:
        """
        Create a new haplotype with a single node entry

        Parameters
        ----------
        node : Variant
            The initializing node for this haplotype
        allele : int
            The allele associated with node

        Returns
        -------
        Haplotype
            The newly created haplotype object containing node and allele
        """
        return cls(((node, allele),))

    def append(self, node: Variant, allele: int) -> Haplotype:
        """
        Append a new node to this haplotype

        Parameters
        ----------
        node : Variant
            The node to add to this haplotype
        allele : int
            The allele associated with this node

        Returns
        -------
        Haplotype
            A new haplotype object extended by the node and its allele
        """
        return Haplotype(self.nodes + ((Variant, allele),))

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
