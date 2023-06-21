#!/usr/bin/env Rscript

# This R script runs the fine-mapping methods SusieR and FINEMAP

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
library(pgenlibr)
library(data.table)


# first, we import the phenotype file into X and Y vectors
args = commandArgs(trailingOnly = TRUE)
gt = args[1]
phen = args[2]
out = args[3]
exclude_causal = as.logical(as.integer(args[4]))

dir.create(out, showWarnings = FALSE)

readPGEN = function(pfile) {
  if (endsWith(pfile, ".pgen")) {
    pfile = substr(pfile, 1, nchar(pfile)-5)
  }
  pvar = pgenlibr::NewPvar(paste0(pfile, '.pvar'))
  pgen = pgenlibr::NewPgen(paste0(pfile, '.pgen'), pvar=pvar)
  n = pgenlibr::GetVariantCt(pgen)
  X = pgenlibr::ReadList(pgen, 1:n, meanimpute=FALSE)
  pgenlibr::ClosePgen(pgen)
  pgenlibr::ClosePvar(pvar)
  X
}

readPSAM = function(pfile) {
  if (endsWith(pfile, ".pgen")) {
    pfile = substr(pfile, 1, nchar(pfile)-5)
  }
  psam = data.table::fread(paste0(pfile, '.psam'), sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE)
  psam[,1]
}

write("reading genotype matrix", stderr())
# import genotype matrices as proper matrices
X = readPGEN(gt)
phen = data.table::fread(phen, sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE)
# the number of samples and the number of variants:
n = nrow(X)
p = ncol(X)
storage.mode(X) = 'double'
# load psam and ensure sample names are the same
stopifnot(readPSAM(gt) == phen[,1])
y = as.matrix(phen[,2])
# what is the column name of the causal variant?
causal_variant = colnames(phen)[2]


# remove the causal variant if requested
if (exclude_causal) {
  X = X[,!(colnames(X) %in% c(causal_variant))]
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
