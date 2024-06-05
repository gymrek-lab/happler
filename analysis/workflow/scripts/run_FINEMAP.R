#!/usr/bin/env Rscript

# This R script runs the fine-mapping method FINEMAP

# param1: The path to a TSV containing the genotype data.
# param2: The path to a TSV containing the phenotype data.
# param3: The path to a directory in which to write output
#         This will be created if it doesn't exist.
# param4: 1 if the causal variant should be removed from the genotype matrix and
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

# first, we import the phenotype file into X and Y vectors
args = commandArgs(trailingOnly = TRUE)
gt = args[1]
phen = args[2]
out = args[3]
exclude_causal = args[4]

dir.create(out, showWarnings = FALSE)

# load functions to help read various file types
source(paste0(thisDir, "/utils.R"))

write("reading genotype matrix", stderr())
# import genotype matrices as proper matrices
X = readPGEN(gt)
phen = readPheno(phen)
# the number of samples and the number of variants:
n = nrow(X)
p = ncol(X)
stopifnot(n > 0)
stopifnot(p > 0)
storage.mode(X) = 'double'
# load psam and ensure sample names are the same
stopifnot(readPSAM(gt) == phen[,1])
y = as.matrix(phen[,2])

# remove the causal variant if requested
causal_variant = NULL
if (exclude_causal != "NULL") {
  stopifnot(readPSAM(exclude_causal) == phen[,1])
  X = cbind(X, readPGEN(exclude_causal))
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
