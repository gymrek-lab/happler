#!/usr/bin/env Rscript

# This R script runs the fine-mapping method FINEMAP

# param1: The path to a TSV containing the genotype data.
# param2: The path to a .glm.linear PLINK2 file containing summary statistics.
# param3: The path to a TSV containing the phenotype data.
# param4: The path to a directory in which to write output
#         This will be created if it doesn't exist.
# param5: 1 if the causal variant should be removed from the genotype matrix and
#         0 otherwise


thisFile <- function() {
  # function to figure out the path to the current script
  # copied from https://stackoverflow.com/a/15373917/16815703
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
thisDir = dirname(thisFile())
write(paste("path to current script is:", thisDir), stderr())


library(abind)
library(data.table)


# first, we import the phenotype file into X and Y vectors
args = commandArgs(trailingOnly = TRUE)
gt = args[1]
phen = args[2]
plink_sumstats = args[3]
out = args[4]
exclude_causal = as.logical(as.integer(args[5]))

dir.create(out, showWarnings = FALSE)


write("reading genotype matrix", stderr())
# import genotype matrices as proper matrices
gt = data.table::fread(gt, sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE)
phen = data.table::fread(phen, sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE)
# create matrices without unecessary columns
# this removes the samples column
X = as.matrix(gt[,-1])
# the number of samples and the number of variants:
n = nrow(X)
p = ncol(X)
storage.mode(X) = 'double'
y = as.matrix(phen[,ncol(phen)])
# what is the column name of the causal variant?
causal_variant = colnames(phen)[2]


# remove the causal variant if requested
if (exclude_causal) {
  X = X[,!(colnames(X) %in% c(causal_variant))]
}


write("computing summary statistics for FINEMAP", stderr())
# compute summary statistics for FINEMAP
mm_regression = function(X, Y, Z=NULL) {
  if (!is.null(Z)) {
      Z = as.matrix(Z)
  }
  reg = lapply(seq_len(ncol(Y)), function (i) simplify2array(susieR:::univariate_regression(X, Y[,i], Z)))
  reg = do.call(abind, c(reg, list(along=0)))
  # return array:
  #   out[1,,] is beta hat (the least-squares estimates of the coefficients)
  #   out[2,,] is se betahat (the standard errors of the beta hats)
  return(aperm(reg, c(3,2,1)))
}
sumstats = mm_regression(X, y)
dat = list(X=X,Y=y)
input = paste0(out,'/sumstats.rds')
write("saving summary statistics to RDS file", stderr())
saveRDS(list(data=dat, sumstats=sumstats), input)

# run FINEMAP
# and set an output path; the results will be written to an RDS file with this basename
output = paste0(out, "/finemap")
args = "--n-causal-snps 1"
commandArgs = function(...) 1
source(paste0(thisDir, "/finemap_1p4.R"))

write("trying to exit with successful error code", stderr())
quit()
