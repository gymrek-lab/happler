from __future__ import annotations
import sys
from pathlib import Path
from typing import TextIO, Generator

import numpy as np
import numpy.typing as npt
from haptools.data import Genotypes

from .tree import Tree
from .variant import Variant


class Haplotype:
    """
    A haplotype within the tree

    Attributes
    ----------
    nodes : tuple[tuple[Variant, int]]
        An ordered collection of pairs, where each pair is a node and its allele
    data : npt.NDArray[np.bool_]
        A np array (with shape n x 2, num_samples x num_chromosomes) denoting the
        presence of this haplotype in each chromosome of each sample
    """

    # TODO: consider using a named tuple?
    nodes: tuple[tuple[Variant, int]]
    data: npt.NDArray[np.bool_]

    def __init__(
        self,
        nodes: tuple[tuple[Variant, int]] = tuple(),
        data: npt.NDArray[np.bool_] = None,
        num_samples: int = None,
    ):
        """
        Initialize an empty haplotype

        Parameters
        ----------
        nodes : tuple[tuple[Variant, int]]
            An ordered collection of pairs, where each pair is a node and its allele
        data : npt.NDArray[np.bool_]
            A np array (with length n x 2, num_samples x num_chromosomes) denoting the
            presence of this haplotype in each chromosome of each sample
        num_samples : int
            The number of samples in this haplotype
        """
        self.nodes = nodes
        if num_samples and data is None:
            self.data = np.ones((num_samples, 2), dtype=np.bool_)
        elif num_samples is None:
            self.data = data
        else:
            raise ValueError(
                "The data and num_samples arguments are mutually exclusive. Provide"
                " either one or the other."
            )

    def __repr__(self):
        return str(self.nodes)

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
            A np array (with length n x 2, num_samples x num_chromosomes) denoting the
            presence of this haplotype in each chromosome of each sample

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
            A np array (with length n x 2, num_samples x num_chromosomes) denoting the
            presence of this haplotype in each chromosome of each sample

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

    def transform(self, genotypes: Genotypes, allele: int) -> npt.NDArray[np.bool_]:
        """
        Transform a genotypes matrix via the current haplotype:

        Each entry in the returned matrix denotes the presence of the current haplotype
        extended by each of the variants in the genotype matrix

        Parameters
        ----------
        genotypes : Genotypes
            The genotypes which to transform using the current haplotype
        allele : int
            The allele (either 0 or 1) of the SNPs we're adding

        Returns
        -------
        npt.NDArray[np.bool_]
            A 3D haplotype matrix similar to the genotype matrix but with haplotypes
            instead of variants in the columns. It will have the same shape except that
            the number of columns (second dimension) will have decreased by the number
            of variants in this haplotype.
        """
        # first, remove any variants that are already in this haplotype using np.delete
        # TODO: consider moving this outside of this function
        gens = np.delete(genotypes.data, self.node_indices, axis=1)
        # add extra axes to match shape of gens
        hap_data = self.data[:, np.newaxis]
        # use np.logical_and to superimpose the current haplotype onto the GT matrix
        return np.logical_and(gens == allele, hap_data)


class Haplotypes:
    """
    A class for processing haplotypes from a file

    Attributes
    ----------
    format : dict
        A dictionary describing the types of lines in the file format as well as
        their format and data types
    data : list[dict]
        A list of dict describing the composition of a series of haplotypes

        Each haplotype dictionary is composed of these items:
            1) id (int): A haplotype ID
            2) tree (int): A tree ID
            3) beta (float): The effect size of the haplotype-phenotype association
            4) pval (float): The p-value of the haplotype-phenotype association
            5) pip (float): A PIP from running SuSiE or some other tool
            6) variants (list[dict]): A list of dictionaries, one for each variant:

        Each variants dictionary is composed of these items:
            1) id (int): A variant ID
            2) hap (int): A haplotype ID
            3) tree (int): A tree ID
            4) allele (bool): The allele for this variant
            5) score (float): The score of this variant within its haplotype
    version : str
        A string denoting the current file format version

    Examples
    --------
    >>> haplotypes = Haplotypes.load('tests/data/simple.haps')
    """

    def __init__(self):
        self.format = {
            "meta": {
                "id": "M",
                "val": ["version"],
                "fmt": ["s"],
            },
            "hap": {
                "id": "H",
                "val": ["id", "tree", "beta", "pval", "pip"],
                "fmt": ["d", "d", ".2f", ".2f", ".2f"],
            },
            "var": {
                "id": "V",
                "val": ["id", "hap", "tree", "allele", "score"],
                "fmt": ["s", "d", "d", "", ".2f"],
            },
        }
        self.version = "0.0.1"
        for val in self.format.keys():
            self.format[val]["str"] = self._create_fmt_str(self.format[val])
        self.data = []

    def _create_fmt_str(self, fmts):
        return (
            fmts["id"]
            + "\t"
            + "\t".join(
                [
                    "{" + val + ":" + fmt + "}"
                    for val, fmt in zip(fmts["val"], fmts["fmt"])
                ]
            )
            + "\n"
        )

    @staticmethod
    def _handle_nan(val, key):
        try:
            if val[key] is not None:
                return val[key]
        except TypeError:
            pass
        return np.nan

    @classmethod
    def from_tree(cls, tree: Tree) -> Haplotypes:
        haps = cls()
        haplotypes = tree.haplotypes()
        haps.data = [
            {
                "id": hap_idx,
                "tree": 0,
                "beta": 0,
                "pval": 0,
                "pip": np.nan,
                "variants": [
                    {
                        "id": node["variant"].id,
                        "hap": hap_idx,
                        "tree": 0,
                        "allele": cls._handle_nan(node, "allele"),
                        "score": cls._handle_nan(node["results"], "pval"),
                    }
                    for node in haplotype
                ],
            }
            for hap_idx, haplotype in enumerate(haplotypes)
        ]
        return haps

    def to_str(self) -> Generator[str, None, None]:
        """
        Create a string representation of this Haplotype

        Yields
        ------
        Generator[str, None, None]
            A list of lines (strings) to include in the output
        """
        yield self.format["meta"]["str"].format(version=self.version)
        for hap in self.data:
            yield self.format["hap"]["str"].format(**hap)
            for var in hap["variants"]:
                yield self.format["var"]["str"].format(**var)

    def __repr__(self):
        return "\n".join(self.to_str())

    def write(self, file: TextIO):
        """
        Write the contents of this Haplotypes object to the file given by fname

        Parameters
        ----------
        file : TextIO
            A file-like object to which this Haplotypes object should be written.
        """
        for line in self.to_str():
            file.write(line)
