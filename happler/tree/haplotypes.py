from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from .tree import Tree
from ..data import Genotypes
from .variant import Variant


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
        gens = np.delete(genotypes.data, self.node_indices, axis=1)
        # add extra axes to match shape of gens
        hap_data = self.data[:, np.newaxis, np.newaxis]
        if gens.shape[1]:
            # use np.logical_and to superimpose the current haplotype onto the GT matrix
            return np.logical_and(gens, hap_data)
        return hap_data


class Haplotypes:
    """
    A class for processing haplotypes from a file

    Attributes
    ----------
    data : npt.NDArray
        An array describing the composition of a series of haplotypes

    Examples
    --------
    >>> haplotypes = Haplotypes.load('tests/data/simple.haps')
    """

    def __init__(self):
        self.data = np.array(
            [],
            dtype=[
                ("hap", np.uint),
                ("tree", np.uint),
                ("variant", "U50"),
                ("allele", np.bool_),
                ("score", np.float64),
            ],
        )

    @classmethod
    def from_tree(cls, tree: Tree) -> Haplotypes:
        haps = cls()
        haplotypes = tree.haplotypes()
        haps.data = np.array(
            [
                (
                    hap_idx,
                    0,
                    node["variant"].id,
                    node["allele"],
                    node["results"]["pval"],
                )
                for hap_idx, haplotype in enumerate(haplotypes)
                for node in haplotype
            ],
            dtype=[
                ("hap", np.uint),
                ("tree", np.uint),
                ("variant", "U50"),
                ("allele", np.bool_),
                ("score", np.float64),
            ],
        )
        return haps
