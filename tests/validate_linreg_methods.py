#!/usr/bin/env python

import numpy as np
from haptools import data
import numpy.typing as npt
from happler.tree.assoc_test import AssocTestSimpleSM as AssocTest

# COMMAND TO TEST AGAINST PLINK2:
# tests/validate_linreg_methods.py > fake.tsv && \
# ~/miniconda3/envs/plink2/bin/plink2 --glm omit-ref allow-no-covars hide-covar \
# --pheno iid-only fake.pheno --pfile fake --no-pheno --out fake &>/dev/null && \
# join -t $'\t' -j1 <(sort -k1,1 fake.tsv) <(cut -f3,9,12 fake.bmi.glm.linear | \
# sort -k1,1) | sort -k3,3g | column -t

sample_size = 100
num_variants = 10
np.random.seed(12345)

gts = data.GenotypesPLINK("fake.pgen")
gts.data = np.random.choice([True, False], size=sample_size*num_variants*2).reshape((sample_size, num_variants, 2))
gts.samples = tuple(f"sample{i}" for i in range(sample_size))
gts.variants = np.array(
	[
		(f"id{i}", "chr0", i, ("A", "T"))
		for i in range(1, num_variants+1)
	],
	dtype=gts.variants.dtype,
)
gts.write()

gts = data.GenotypesPLINK(gts.fname)
gts.read()

def standardize(X: npt.NDArray[np.uint8]) -> npt.NDArray[np.float64]:
    """
    Standardize the genotypes so they have mean 0 and variance 1

    Parameters
    ----------
    X : npt.NDArray[np.float64]
        The genotypes, with shape n x p. There are only two dimensions.
        Each column is a haplotype and each row is a sample.

    Returns
    -------
    npt.NDArray[np.float64]
        An array with the same shape as X but standardized properly
    """
    std = np.std(X, axis=0)
    standardized = (X - np.mean(X, axis=0)) / std
    # # for variants where the stdev is 0, just set all values to 0 instead of nan
    # zero_elements = std == 0
    # standardized[:, zero_elements] = np.zeros(
    #     (X.shape[0], np.sum(zero_elements))
    # )
    return standardized

pts = data.Phenotypes("fake.pheno")
# pts.data = np.random.normal(size=gts.data.shape[0]) * 0.4
pts.data = gts.data[:, 0].sum(axis=1)*0.005 + np.random.normal(scale=0.6, size=gts.data.shape[0]) + 2
pts.data = pts.data[:, np.newaxis]
pts.samples = tuple(gts.samples)
pts.names = ("bmi",)
pts.write()

pts = data.Phenotypes(pts.fname)
pts.read()
pts.subset(samples=gts.samples, inplace=True)
pts.standardize()

tester = AssocTest()
results = tester.run(gts.data[:, :, :2].sum(axis=2), pts.data[:,0].flatten()).data
# for i in range(num_variants):
# for i in [0]:
# 	plt.scatter(gts.data[:, i].sum(axis=1), pts.data)
# 	plt.savefig(f"id{i+1}.png")
# 	plt.clf()
print("ID\tBETA\tP")
for ID, beta, pval in zip(gts.variants["id"], results["beta"], results["pval"]):
	print(f"{ID}\t{beta}\t{pval}")
