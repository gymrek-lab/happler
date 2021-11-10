import numpy as np
from cyvcf2 import VCF, Variant
from pathlib import Path
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
        self.loaded = False
        self.fname = fname
        self.data = np.array([])
        super().__init__()

    @abstractmethod
    def load(self):
        """
        Load the file contents into self.data
        """
        pass


class Genotypes(Data):
    """
    A class for reading genotypes

    Attributes
    ----------
    data : np.array
        The genotypes in an n x p x 2 array
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
        Read genotypes from a VCF into a numpy matrix
        """
        if self.loaded:
            raise AssertionError(
                "The data has already been loaded."
            )
        # load all info into memory
        vcf = VCF(str(self.fname))
        self.samples = vcf.samples
        variants = list(vcf)
        # save meta information about each variant
        self.variants = np.array([
            (variant.ID, variant.CHROM, variant.POS, variant.aaf)
            for variant in variants
        ], dtype=[
            ('id', 'U50'), ('chrom', 'U10'), ('pos', np.uint), ('aaf', np.float64)
        ])
        # extract the genotypes to a np matrix of size n x p x 3
        # the last dimension has three items:
        # 1) presence of REF in strand one
        # 2) presence of REF in strand two
        # 3) whether the genotype is phased
        self.data = np.array([
            variant.genotypes for variant in variants
        ], dtype=np.bool_)
        self.loaded = True

    def check_phase(self):
        """
        Check that the genotypes are phased then remove the phasing info from the data
        """
        if self.data.shape[2] < 3:
            raise AssertionError(
                "Phase information has already been removed from the data"
            )
        # check: are there any variants that are heterozygous and unphased?
        unphased = (self.data[:,:,0] ^ self.data[:,:,1]) & (self.data[:,:,2])
        if np.any(unphased):
            variant_idx, samp_idx = np.nonzero(unphased)
            raise ValueError(
                "Variant with ID {} at POS {}:{} is unphased for sample {}".format(
                    *self.variants[variant_idx], vcf.samples[samp_idx]
                )
            )
        # remove the last dimension that contains the phase info
        self.data = self.data[:,:,:2]

    def to_MAC(self):
        """
        convert an ALT count GT matrix into a matrix of minor allele counts
        """
        if self.mac_converted:
            raise AssertionError(
                "The matrix already counts instances of the minor allele rather than"
                "the ALT allele."
            )
        need_conversion = self.variants['aaf'] > 0.5
        # flip the strands on the variants that have an alternate allele frequency
        # above 0.5
        self.data = self.data[need_conversion, :, ::-1]
        self.variants['aaf'][need_conversion] = 1 - self.variants['aaf'][need_conversion]
        dtype = list(self.variants.dtype.names)
        dtype[-1] = 'maf'
        self.variants.astype(dtype)
