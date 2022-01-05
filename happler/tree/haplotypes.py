from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

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
    data : npt.NDArray[np.bool_]
        A np array (with shape n x 1, the number of samples) denoting the presence
        of this haplotype in each sample
    """

    # TODO: consider using a named tuple?
    nodes: tuple[tuple[Variant, int]]
    data: npt.NDArray[np.bool_]

    @classmethod
    def from_node(
        cls, node: Variant, allele: int, variant_genotypes: npt.NDArray[np.bool_]
    ) -> Haplotype:
        """
        Create a new haplotype with a single node entry

        Parameters
        ----------
        node : Variant
            The initializing node for this haplotype
        allele : int
            The allele associated with node
        variant_genotypes : npt.NDArray[np.bool_]
            A np array (with shape n x 1, the number of samples) denoting the presence of
            this genotype in each sample

        Returns
        -------
        Haplotype
            The newly created haplotype object containing node and allele
        """
        return cls(((node, allele),), variant_genotypes)

    def append(
        self, node: Variant, allele: int, variant_genotypes: npt.NDArray[np.bool_]
    ) -> Haplotype:
        """
        Append a new node (variant) to this haplotype

        Parameters
        ----------
        node : Variant
            The node to add to this haplotype
        allele : int
            The allele associated with this node
        variant_genotypes : npt.NDArray[np.bool_]
            A np array (with length n x 1, the number of samples) denoting the presence of
            this genotype in each sample

        Returns
        -------
        Haplotype
            A new haplotype object extended by the node and its allele
        """
        new_haplotype = self.nodes + ((node, allele),)
        new_haplotype_values = self.data & variant_genotypes
        return Haplotype(new_haplotype, new_haplotype_values)

    # TODO: use @cached_property from python3.8, instead
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

    def transform(self, genotypes: Genotypes) -> npt.NDArray[np.bool_]:
        """
        Transform a genotypes matrix via the current haplotype:

        Each entry in the returned matrix denotes the presence of the current haplotype
        extended by each of the variants in the genotype matrix

        Parameters
        ----------
        genotypes : Genotypes
            The genotypes which to transform using the current haplotype

        Returns
        -------
        npt.NDArray[np.bool_]
            A haplotype matrix similar to the genotype matrix but with haplotypes
            instead of variants in the columns. It will have the same shape except that
            the number of columns (second dimension) will have decreased by one.
        """
        # first, remove any variants that are already in this haplotype using np.delete
        # then, use a np.logical_and to impose the current haplotype onto the GT matrix
        return np.logical_and(
            np.delete(genotypes.data, self.node_indices, axis=1),
            self.data[:, np.newaxis, np.newaxis],
        )


class Haplotypes:
    """
    A class for processing haplotypes from a file

    Attributes
    ----------
    genotypes : Genotypes
        The original genotypes from which these Haplotypes were derived
    data : np.array
        The haplotypes in an n (samples) x p (variants) x 2 (strands) array
    haplotypes : list[Haplotype]
        Haplotype-level meta information

    Examples
    --------
    >>> haplotypes = Haplotypes.load('tests/data/simple.vcf')
    """

    def __init__(self, genotypes: Genotypes):
        self.genotypes = genotypes
        # initialize data to a matrix composed of no haplotypes
        self.data = np.empty((len(self.genotypes.samples), 0), dtype=np.bool_)
        self.haplotypes = []

    def add_variant(self, variant: Variant, allele: int, hap_idxs: list[int] = None):
        # do we already have any haplotypes?
        variant_gts = self.genotypes.data[:, variant.idx, allele][:, np.newaxis]
        if self.haplotypes:
            if hap_idxs is None:
                hap_idxs = list(range(len(self.haplotypes)))
            for hap_idx in hap_idxs:
                self.haplotypes[hap_idx] = self.haplotypes[hap_idx].append(
                    variant, allele, variant_gts
                )
            self.data[:, hap_idxs] = np.logical_and(self.data[:, hap_idxs], variant_gts)
        else:
            self.haplotypes = [Haplotype.from_node(variant, allele, variant_gts)]
            self.data = variant_gts

    # WHAT DO WE NEED THIS CLASS TO DO?
    # 1) create a haplotype matrix where each (SNP, allele) pair is added to an existing haplotype
    # - at the root of the tree, this is just the same as add_variant()
    # 2) compose a vector of haplotype dosages from the haplotype matrix
