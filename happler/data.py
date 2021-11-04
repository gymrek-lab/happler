import numpy as np
from cyvcf2 import VCF
from pathlib import Path

from IPython import embed


class Data:
    """
    Class for accessing read-only data
    """

    def __init__(self, gens: Path, phens: Path):
        """
        Read genotype and phenotype files

        :param gens: The path to the genotypes
        :param phens: The path to the phenotypes
        """
        self.gens = self.read_genotypes(gens)
        self.phens = None # TODO

    def read_genotypes(self, gens: Path):
        """
        Read genotypes from a VCF into a numpy matrix

        :param gens: The path to the VCF file
        :return: A genotype matrix; rows are samples and columns are SNP genotypes
        :rtype: np.array
        """
        vcf = VCF(str(gens))
        embed()
        num_samps = len(vcf.samples)
        # extract the start, REF, and ALT alleles of each variant
        gts = np.array([
            (variant.start, *tuple(sum_gts(variant.genotypes)))
            for variant in VCF
        ], dtype='|i8'+(', U1'*num_samps))
        return gts

    def sum_gts(self, gts: list):
        """
        Sum the genotypes across a list of genotypes.
        For example, [0, 1] --> 1

        :param gts: The genotypes
        :return: A list of summed genotypes
        :type: list
        """
        for gt in gts:
            yield sum(gt[:2])
