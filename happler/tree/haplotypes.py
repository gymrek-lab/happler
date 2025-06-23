from __future__ import annotations
import sys
from pathlib import Path
from logging import Logger
from typing import TextIO, Generator
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
from haptools.data import (
    Extra,
    Genotypes,
    GenotypesVCF,
    Variant as VariantBase,
    Haplotype as HaplotypeBase,
    Haplotypes as HaplotypesBase,
)

from .tree import Tree
from .variant import Variant


class Haplotype:
    """
    A haplotype within the tree

    Attributes
    ----------
    nodes : tuple[tuple[Variant, int]]
        An ordered collection of pairs, where each pair is a node and its allele
    data : npt.NDArray[bool]
        A np array (with shape n x 2, num_samples x num_chromosomes) denoting the
        presence of this haplotype in each chromosome of each sample
    """

    # TODO: consider using a named tuple?
    nodes: tuple[tuple[Variant, int]]
    data: npt.NDArray[bool]

    def __init__(
        self,
        nodes: tuple[tuple[Variant, int]] = tuple(),
        data: npt.NDArray[bool] = None,
        num_samples: int = None,
    ):
        """
        Initialize an empty haplotype

        Parameters
        ----------
        nodes : tuple[tuple[Variant, int]]
            An ordered collection of pairs, where each pair is a node and its allele
        data : npt.NDArray[bool]
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
        cls, node: Variant, allele: int, variant_genotypes: npt.NDArray[bool]
    ) -> Haplotype:
        """
        Create a new haplotype with a single node entry

        Parameters
        ----------
        node : Variant
            The initializing node for this haplotype
        allele : int
            The allele associated with node
        variant_genotypes : npt.NDArray[bool]
            A np array (with length n x 2, num_samples x num_chromosomes) denoting the
            presence of this haplotype in each chromosome of each sample

        Returns
        -------
        Haplotype
            The newly created haplotype object containing node and allele
        """
        return cls(((node, allele),), variant_genotypes)

    @classmethod
    def from_haptools_haplotype(
        cls, haplotype: HaplotypeBase, variant_genotypes: GenotypesVCF
    ) -> Haplotype:
        """
        Create a new haplotype from a haptools Haplotype and a GenotypesVCF object
        """
        variants = {
            vr.id: vr.allele for vr in haplotype.variants
        }
        variant_genotypes.index(variants = True)
        gts = variant_genotypes.subset(variants=tuple(variants.keys()))
        nodes = tuple(
            (
                Variant.from_np(variant, variant_genotypes._var_idx[variant["id"]]),
                list(variant["alleles"]).index(allele)
            )
            for variant, allele in zip(gts.variants, variants.values())
        )
        return cls(nodes, haplotype.transform(gts))

    def append(
        self, node: Variant, allele: int, variant_genotypes: npt.NDArray[bool]
    ) -> Haplotype:
        """
        Append a new node (variant) to this haplotype

        Parameters
        ----------
        node : Variant
            The node to add to this haplotype
        allele : int
            The allele associated with this node
        variant_genotypes : npt.NDArray[bool]
            A np array (with length n x 2, num_samples x num_chromosomes) denoting the
            presence of the new allele in each chromosome of each sample

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

    def transform(
        self,
        genotypes: Genotypes,
        allele: int,
        idxs: tuple[int] = None,
    ) -> npt.NDArray[bool]:
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
        idxs : tuple[int], optional
            If specified, we will only output haplotypes for the variants at these
            indices. Otherwise, we'll output all of them.

        Returns
        -------
        npt.NDArray[bool]
            A 3D haplotype matrix similar to the genotype matrix but with haplotypes
            instead of variants in the columns. It will have the same shape except that
            the number of columns (second dimension) will have decreased by the number
            of variants in this haplotype.
        """
        # first, remove any variants that are already in this haplotype using np.delete
        # TODO: consider moving this outside of this function
        gens = np.delete(genotypes.data, self.node_indices, axis=1)
        # how does the deletion change the desired indices?
        if idxs is not None:
            idxs -= np.sum(np.array(self.node_indices)[:, np.newaxis] < idxs, axis=0)
        else:
            # alias for all of the indices
            idxs = np.s_[:]
        # add extra axes to match shape of gens
        hap_data = self.data[:, np.newaxis]
        # use np.logical_and to superimpose the current haplotype onto the GT matrix
        return np.logical_and(gens[:, idxs] == allele, hap_data)


@dataclass
class HapplerVariant(VariantBase):
    """
    A variant allele with sufficient fields for happler
    Properties and functions are shared with the base Variant object, "VariantBase"
    """

    score: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(Extra("score", ".2f", "BIC assigned to this variant"),),
    )


@dataclass
class HapplerHaplotype(HaplotypeBase):
    """
    A haplotype with sufficient fields for happler
    Properties and functions are shared with the base Haplotype object, "HaplotypeBase"
    """

    beta: float
    pval: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(
            Extra("beta", ".2f", "Effect size in linear model"),
            Extra("pval", ".2f", "-log(pval) in linear model"),
        ),
    )


class Haplotypes(HaplotypesBase):
    """
    A class for processing haplotypes from a file

    Attributes
    ----------
    fname: Path | str
        The path to the file containing the data
    data: dict[str, Haplotype]
        A dict of Haplotype objects keyed by their IDs
    types: dict
        A dict of class names keyed by the symbol denoting their line type
        Ex: {'H': Haplotype, 'V': Variant}
    version: str
        A string denoting the current file format version
    log: Logger
        A logging instance for recording debug statements.
    """

    @staticmethod
    def _handle_nan(val, key):
        try:
            if val[key] is not None:
                return val[key]
        except TypeError:
            pass
        return np.nan

    @classmethod
    def from_tree(
        cls, fname: Path | str, tree: Tree, gts: GenotypesVCF, log: Logger = None
    ) -> Haplotypes:
        """
        Create a Haplotypes object from a Tree object and a Genotypes object

        Parameters
        ----------
        fname : Path | str
            The fname parameter for the Haplotypes object
        tree : Tree
            The Tree object containing the haplotypes to encode within a Haplotypes obj
        gts : GenotypesVCF
            The genotypes from which the tree was constructed
        log : Logger, optional
            The log parameter for the Haplotypes object

        Returns
        -------
        Haplotypes
            The completed Haplotypes object
        """
        haps = cls(
            fname=fname, haplotype=HapplerHaplotype, variant=HapplerVariant, log=log
        )
        haps.data = {}
        for hap_idx, haplotype in enumerate(tree.haplotypes()):
            hap_id = "H" + str(hap_idx)
            results = haplotype[-1]["results"]
            haps.data[hap_id] = HapplerHaplotype(
                chrom=gts.variants[haplotype[0]["variant"].idx]["chrom"],
                start=0,  # this is filled out later
                end=0,  # this is filled out later
                id=hap_id,
                beta=results["beta"],
                pval=-np.log10(results["pval"]),
            )
            alleles = {
                node["variant"].idx: gts.variants[node["variant"].idx]["alleles"][
                    cls._handle_nan(node, "allele")
                ]
                for node in haplotype
            }
            haps.data[hap_id].variants = tuple(
                HapplerVariant(
                    start=node["variant"].pos,
                    end=node["variant"].pos + len(alleles[node["variant"].idx]),
                    id=node["variant"].id,
                    allele=alleles[node["variant"].idx],
                    score=cls._handle_nan(node["results"], "bic"),
                )
                for node in haplotype
            )
            haps.data[hap_id].start = min(n.start for n in haps.data[hap_id].variants)
            haps.data[hap_id].end = max(n.end for n in haps.data[hap_id].variants)
        return haps
