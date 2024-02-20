#!/usr/bin/env Rscript

# This R script creates PEER factors for expression data

# param1: The path to a CSV containing a samples by genes matrix of phenotypes (ex: EUR_converted_expr_hg38.csv)
# param2: The path to a CSV to which the output should be written

library(peer)

# first, we import the necessary files
args = commandArgs(trailingOnly = TRUE)
phen = args[1]
out = args[2] # ex: EUR_normalized_filtered_peer_factors_hg38.csv

# Input should be Samples x Genes (rows x cols)
write("reading phenotype matrix", stderr())
expr = read.csv(phen, header=FALSE)
write(dim(expr), stderr())
model = PEER()
write("setting pheno model", stderr())
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setNk(model,36) # Set number of Peer factors
write("running training process", stderr())
PEER_update(model) # Run training process
# residuals = PEER_getResiduals(model) # Calc Residuals
# dim(residuals)
write("retrieving PEER factors", stderr())
factors = PEER_getX(model)
dim(factors)
# Write factors
write("writing PEER factors", stderr())
write.table(factors, file=out, sep=" ", col.names=TRUE, row.names=TRUE)
