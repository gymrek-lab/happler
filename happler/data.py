import numpy as np
from csv import reader
from pathlib import Path
from cyvcf2 import VCF, Variant
from abc import ABC, abstractmethod


class Data(ABC):
    """
    Abstract class for accessing read-only data files

    Attributes
    ----------
    data : np.array
        The contents of the data file, once loaded
    fname : Path
        The path to the read-only file containing the data
    """

    def __init__(self, fname: Path):
        self.fname = fname
        self.data = None
        super().__init__()

    @abstractmethod
    def load(self):
        """
        Load the file contents into self.data
        """
        if self.data is not None:
            raise AssertionError("The data has already been loaded.")


class Genotypes(Data):
    """
    A class for reading genotypes

    Attributes
    ----------
    data : np.array
        The genotypes in an n (samples) x p (variants) x 2 (strands) array
    samples : list
        Sample-level meta information
    variants : list
        Variant-level meta information
    """

    def __init__(self, fname: Path):
        super().__init__(fname)
        self.samples = []
        self.variants = np.array([])

    def load(self):
        """
        Read genotypes from a VCF into a numpy matrix stored in `self.data`
        """
        super().load()
        # load all info into memory
        vcf = VCF(str(self.fname))
        self.samples = vcf.samples
        variants = list(vcf)
        # save meta information about each variant
        self.variants = np.array(
            [
                (variant.ID, variant.CHROM, variant.POS, variant.aaf)
                for variant in variants
            ],
            dtype=[
                ("id", "U50"),
                ("chrom", "U10"),
                ("pos", np.uint),
                ("aaf", np.float64),
            ],
        )
        # extract the genotypes to a np matrix of size n x p x 3
        # the last dimension has three items:
        # 1) presence of REF in strand one
        # 2) presence of REF in strand two
        # 3) whether the genotype is phased
        self.data = np.array(
            [variant.genotypes for variant in variants], dtype=np.bool_
        )
        # transpose the GT matrix so that samples are rows and variants are columns
        self.data = self.data.transpose((1, 0, 2))

    def check_phase(self):
        """
        Check that the genotypes are phased then remove the phasing info from the data

        Raises
        ------
        AssertionError
            If the phase information has already been checked and removed from the data
        ValueError
            If any heterozgyous genotpyes are unphased
        """
        if self.data.shape[2] < 3:
            raise AssertionError(
                "Phase information has already been removed from the data"
            )
        # check: are there any variants that are heterozygous and unphased?
        unphased = (self.data[:, :, 0] ^ self.data[:, :, 1]) & (~self.data[:, :, 2])
        if np.any(unphased):
            samp_idx, variant_idx = np.nonzero(unphased)
            raise ValueError(
                "Variant with ID {} at POS {}:{} is unphased for sample {}".format(
                    *tuple(self.variants[variant_idx[0]])[:3], self.samples[samp_idx[0]]
                )
            )
        # remove the last dimension that contains the phase info
        self.data = self.data[:, :, :2]

    def to_MAC(self):
        """
        Convert the ALT count GT matrix into a matrix of minor allele counts

        Raises
        ------
        AssertionError
            If the matrix has already been converted
        """
        if self.variants.dtype.names[3] == "maf":
            raise AssertionError(
                "The matrix already counts instances of the minor allele rather than"
                "the ALT allele."
            )
        need_conversion = self.variants["aaf"] > 0.5
        # flip the strands on the variants that have an alternate allele frequency
        # above 0.5
        self.data[:, need_conversion, :2] = self.data[:, need_conversion, 1::-1]
        # also encode an MAF instead of an AAF in self.variants
        self.variants["aaf"][need_conversion] = (
            1 - self.variants["aaf"][need_conversion]
        )
        # replace 'aaf' with 'maf' in the matrix
        self.variants.dtype.names = [
            (x, "maf")[x == "aaf"] for x in self.variants.dtype.names
        ]
