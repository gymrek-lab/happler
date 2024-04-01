library(pgenlibr)
library(data.table)

readPVAR_base = function(pfile) {
  if (endsWith(pfile, ".pgen")) {
    pfile = substr(pfile, 1, nchar(pfile)-5)
  }
  pvar = data.table::fread(paste0(pfile, '.pvar'), sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE, skip="#CHROM\tPOS")
  colnames(pvar)[1] = "CHROM"
  pvar$POS = as.integer(pvar$POS)
  pvar
}

readPVAR = function(pfile, region=NULL) {
  pvar = readPVAR_base(pfile)
  if (!is.null(region)) {
    chrom = strsplit(region, split=":")[[1]][1]
    start_end = strsplit(strsplit(region, split=":")[[1]][2], split="-")[[1]]
    matching = pvar$CHROM == chrom & pvar$POS > start_end[1] & pvar$POS < start_end[2]
    pvar = pvar[matching,]
  }
  pos = as.vector(pvar$POS)
  names(pos) = pvar$ID
  pos
}

readPVAR_region = function(pfile, region) {
  chrom = strsplit(region, split=":")[[1]][1]
  start_end = strsplit(strsplit(region, split=":")[[1]][2], split="-")[[1]]
  pvar = readPVAR_base(pfile)
  which(pvar$CHROM == chrom & pvar$POS > start_end[1] & pvar$POS < start_end[2])
}

readPGEN = function(pfile, region=NULL) {
  if (endsWith(pfile, ".pgen")) {
    pfile = substr(pfile, 1, nchar(pfile)-5)
  }
  pvar = pgenlibr::NewPvar(paste0(pfile, '.pvar'))
  pgen = pgenlibr::NewPgen(paste0(pfile, '.pgen'), pvar=pvar)
  if (is.null(region)) {
    n = pgenlibr::GetVariantCt(pgen)
    var_subset = 1:n
  } else {
    var_subset = readPVAR_region(pfile, region)
  }
  X = pgenlibr::ReadList(pgen, var_subset, meanimpute=FALSE)
  pgenlibr::ClosePgen(pgen)
  pgenlibr::ClosePvar(pvar)
  colnames(X) = names(readPVAR(pfile, region=region))
  X
}

readPSAM = function(pfile) {
  if (endsWith(pfile, ".pgen")) {
    pfile = substr(pfile, 1, nchar(pfile)-5)
  }
  psam = data.table::fread(paste0(pfile, '.psam'), sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE)
  psam[,1]
}

readPheno = function(pheno) {
  data.table::fread(pheno, sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE)
}

save_curr_env = function(env = parent.frame(), filePath = "myenv.RData") {
  # This function is useful for debugging failed R scripts
  # When this function is called from within a script, it will save all current
  # environment variables into the specified RData file.
  # Parameters:
  #   env: The environment to save. By default, it saves the environment from which
  #        the function was called (parent.frame()).
  #   filePath: The path of the file where the environment will be saved.
  #              By default, it saves to "myenv.RData" in the current working directory.
  #
  # Within an R session on the command line, you can then execute load("myenv.RData")
  # or load(filePath) to inspect the environment.

  save(list = ls(all.names = TRUE, envir = env), file = filePath, envir = env)
}
