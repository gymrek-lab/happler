#!/usr/bin/env python

# COMMAND TO TEST AGAINST PLINK2:

# tests/validate_linreg_methods.py > fake.tsv && \
# ~/miniconda3/envs/plink2/bin/plink2 --glm omit-ref hide-covar --covar iid-only fake.covar \
# --pheno iid-only fake.pheno --pfile fake --threads 1 --no-pheno --out fake &>/dev/null && \
# diff -sy -W 67 <(sort -k1,1 fake.tsv | column -t) <(cut -f3,9,12 fake.bmi.glm.linear | sort -k1,1 | column -t)


import numpy as np
from haptools import data
import numpy.typing as npt
from happler.tree.assoc_test import AssocTestSimpleCovariates as AssocTest

sample_size = 100
num_variants = 5
np.random.seed(12345)

gts = data.GenotypesPLINK("fake.pgen")
gts.data = np.random.choice([True, False], size=sample_size * num_variants * 2).reshape(
    (sample_size, num_variants, 2)
)
gts.samples = tuple(f"sample{i}" for i in range(sample_size))
gts.variants = np.array(
    [(f"id{i}", "chr0", i, ("A", "T")) for i in range(1, num_variants + 1)],
    dtype=gts.variants.dtype,
)
gts.write()

gts = data.GenotypesPLINK(gts.fname)
gts.read()

pts = data.Phenotypes("fake.pheno")
# pts.data = np.random.normal(size=gts.data.shape[0]) * 0.4
pts.data = (
    gts.data[:, 0].sum(axis=1) * 0.4
    + np.random.normal(scale=0.6, size=gts.data.shape[0])
    + 2
)
pts.data = pts.data[:, np.newaxis]
pts.samples = tuple(gts.samples)
pts.names = ("bmi",)
pts.write()

pts = data.Phenotypes(pts.fname)
pts.read()
pts.subset(samples=gts.samples, inplace=True)
pts.standardize()

cvs = data.Covariates("fake.covar")
cvs.data = np.random.normal(size=gts.data.shape[0]) * 0.01 + 2
cvs.data = cvs.data[:, np.newaxis]
cvs.samples = tuple(gts.samples)
cvs.names = ("sex",)
cvs.write()

cvs = data.Covariates(cvs.fname)
cvs.read()
cvs.subset(samples=gts.samples, inplace=True)
cvs.standardize()

tester = AssocTest(covars=cvs.data)
results = tester.run(gts.data[:, :, :2].sum(axis=2), pts.data[:, 0].flatten()).data
# for i in range(num_variants):
# for i in [0]:
# 	plt.scatter(gts.data[:, i].sum(axis=1), pts.data)
# 	plt.savefig(f"id{i+1}.png")
# 	plt.clf()
print("ID\tBETA\tP")
for ID, beta, pval in zip(gts.variants["id"], results["beta"], results["pval"]):
    print(f"{ID}\t{beta:.6g}\t{pval:.6g}")
