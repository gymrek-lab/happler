from __future__ import annotations
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

    def __repr__(self):
        return str(self.fname)

    @classmethod
    @abstractmethod
    def load(cls: Data, fname: Path):
        """
        Read the file contents and perform any recommended pre-processing

        Parameters
        ----------
        fname : Path
            See documentation for :py:attr:`~.Data.fname`
        """
        pass

    @abstractmethod
    def read(self):
        """
        Read the raw file contents into the class properties
        """
        if self.data is not None:
            raise AssertionError("The data has already been loaded.")


class Genotypes(Data):
    """
    A class for processing genotypes from a file

    Attributes
    ----------
    data : np.array
        The genotypes in an n (samples) x p (variants) x 2 (strands) array
    samples : tuple
        The names of each of the n samples
    variants : list
        Variant-level meta information:
            1. ID
            2. CHROM
            3. POS
            4. AAF: allele freq of alternate allele (or MAF if to_MAC() is called)

    Examples
    --------
    >>> genotypes = Genotypes.load('tests/data/simple.vcf')
    """

    def __init__(self, fname: Path):
        super().__init__(fname)
        self.samples = tuple()
        self.variants = np.array([])

    @classmethod
    def load(cls: Genotypes, fname: Path) -> Genotypes:
        """
        Load genotypes from a VCF file

        Read the file contents, check the genotype phase, and create the MAC matrix

        Parameters
        ----------
        fname
            See documentation for :py:attr:`~.Data.fname`

        Returns
        -------
        genotypes
            A Genotypes object with the data loaded into its properties
        """
        genotypes = cls(fname)
        genotypes.read()
        genotypes.check_biallelic()
        genotypes.check_phase()
        # genotypes.to_MAC()
        return genotypes

    def read(self):
        """
        Read genotypes from a VCF into a numpy matrix stored in :py:attr:`~.Genotypes.data`
        """
        super().read()
        # load all info into memory
        vcf = VCF(str(self.fname))
        self.samples = tuple(vcf.samples)
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
            [variant.genotypes for variant in variants], dtype=np.uint8
        )
        # transpose the GT matrix so that samples are rows and variants are columns
        self.data = self.data.transpose((1, 0, 2))

    def check_biallelic(self):
        """
        Check that each genotype is composed of only two alleles

        This function modifies the dtype of :py:attr:`~.Genotypes.data` from uint8 to bool

        Raises
        ------
        AssertionError
            If the number of alleles has already been checked and the dtype has been
            converted to bool
        ValueError
            If any of the genotypes have more than two alleles
        """
        if self.data.dtype == np.bool_:
            raise AssertionError("All genotypes are already biallelic")
        # check: are there any variants that have genotype values above 1?
        # A genotype value above 1 would imply the variant has more than one ALT allele
        multiallelic = np.any(self.data[:, :, :2] > 1, axis=2)
        if np.any(multiallelic):
            samp_idx, variant_idx = np.nonzero(multiallelic)
            raise ValueError(
                "Variant with ID {} at POS {}:{} is multiallelic for sample {}".format(
                    *tuple(self.variants[variant_idx[0]])[:3], self.samples[samp_idx[0]]
                )
            )
        self.data = self.data.astype(np.bool_)

    def check_phase(self):
        """
        Check that the genotypes are phased then remove the phasing info from the data

        This function modifies :py:attr:`~.Genotypes.data` in-place

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

        This function modifies :py:attr:`~.Genotypes.data` in-place

        It also changes the 'aaf' record in :py:attr:`~.Genotypes.variants` to 'maf'

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
        # flip the count on the variants that have an alternate allele frequency
        # above 0.5
        self.data[:, need_conversion, :2] = ~self.data[:, need_conversion, :2]
        # also encode an MAF instead of an AAF in self.variants
        self.variants["aaf"][need_conversion] = (
            1 - self.variants["aaf"][need_conversion]
        )
        # replace 'aaf' with 'maf' in the matrix
        self.variants.dtype.names = [
            (x, "maf")[x == "aaf"] for x in self.variants.dtype.names
        ]


class Phenotypes(Data):
    """
    A class for processing phenotypes from a file

    Attributes
    ----------
    data : np.array
        The phenotypes in an n (samples) x 1 (phenotype value) array
    samples : tuple
        The names of each of the n samples

    Examples
    --------
    >>> phenotypes = Phenotypes.load('tests/data/simple.tsv')
    """

    def __init__(self, fname: Path):
        super().__init__(fname)
        self.samples = tuple()

    @classmethod
    def load(cls: Phenotypes, fname: Path) -> Phenotypes:
        """
        Load phenotypes from a TSV file

        Read the file contents and standardize the phenotypes

        Parameters
        ----------
        fname
            See documentation for :py:attr:`~.Data.fname`

        Returns
        -------
        phenotypes
            A Phenotypes object with the data loaded into its properties
        """
        phenotypes = cls(fname)
        phenotypes.read()
        phenotypes.standardize()
        return phenotypes

    def read(self):
        """
        Read phenotypes from a TSV file into a numpy matrix stored in :py:attr:`~.Penotypes.data`

        Raises
        ------
        AssertionError
            If the provided file doesn't follow the expected format
        """
        super().read()
        # load all info into memory
        with open(self.fname) as phens:
            phen_text = list(reader(phens, delimiter="\t"))
        # there should only be two columns
        assert len(phen_text[0]) == 2, "The phenotype TSV should only have two columns."
        # the second column should be castable to a float
        try:
            float(phen_text[0][1])
        except:
            raise AssertionError("The second column of the TSV file must numeric.")
        # fill out the samples and data properties
        self.samples, self.data = zip(*phen_text)
        # coerce strings to floats
        self.data = np.array(self.data, dtype="float64")

    def standardize(self):
        """
        Standardize phenotypes so they have a mean of 0 and a stdev of 1

        This function modifies :py:attr:`~.Genotypes.data` in-place
        """
        self.data = (self.data - np.mean(self.data)) / np.std(self.data)
