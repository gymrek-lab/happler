#!/usr/bin/env Rscript

# This R script runs the fine-mapping method SusieR

# param1: The path to a TSV containing the genotype data.
# param2: The path to a TSV containing the phenotype data.
# param3: The path to a directory in which to write output
#         This will be created if it doesn't exist.
# param4: 1 if the causal variant should be removed from the genotype matrix and
#         0 otherwise
# param5: The finemapping region (ex: 19:45401409-46401409)


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
region = args[5]

dir.create(out, showWarnings = FALSE)

# load functions to help read various file types
source(paste0(thisDir, "/utils.R"))

write("reading genotype matrix and phenotypes", stderr())
samples = readPSAM(gt, samples=readPheno(phen)[,1])
phen = readPheno(phen, samples=samples[,1])
# load psam and ensure sample names are the same
stopifnot(samples[,1] == phen[,1])
y = as.matrix(phen[,2])
# import genotype matrices as proper matrices
X = readPGEN(gt, region=region, samples=samples[,2])
# the number of samples and the number of variants:
n = nrow(X)
p = ncol(X)
stopifnot(n > 0)
stopifnot(p > 0)
storage.mode(X) = 'double'

# remove the causal variant if requested
causal_variant = NULL
if (exclude_causal != "NULL") {
  exclude_causal_samples = readPSAM(exclude_causal, samples=phen[,1])
  stopifnot(exclude_causal_samples[,1] == phen[,1])
  X = cbind(X, readPGEN(exclude_causal, samples=exclude_causal_samples[,2]))
}

# run SuSiE
# write the output to an RDS file
write("executing SuSiE", stderr())
fitted = susieR::susie(X, y, L=1)

# when writing the output, also include information about which variant is causal
# and whether it was included in the simulation
write("writing SuSiE output", stderr())
saveRDS(
  list(causal_var=causal_variant, causal_excluded=exclude_causal, fitted=fitted),
  paste0(out, '/susie.rds')
)

write("trying to exit with successful error code", stderr())
quit()
