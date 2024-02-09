library(pgenlibr)
library(data.table)

readPVAR = function(pfile) {
  if (endsWith(pfile, ".pgen")) {
    pfile = substr(pfile, 1, nchar(pfile)-5)
  }
  pvar = data.table::fread(paste0(pfile, '.pvar'), sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE, skip="#CHROM\tPOS")
  pos = as.integer(as.vector(pvar$POS))
  names(pos) = pvar$ID
  pos
}

readPVAR_region = function(pfile, region) {
  chrom = strsplit(region, split=":")[1]
  start_end = strsplit(strsplit(region, split=":")[2], split="-")
  if (endsWith(pfile, ".pgen")) {
    pfile = substr(pfile, 1, nchar(pfile)-5)
  }
  pvar = data.table::fread(paste0(pfile, '.pvar'), sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE, skip="#CHROM\tPOS")
  which(pvar$CHROM == chrom && pvar$POS < start_end[1] && pvar$POS > start_end[2])
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
  colnames(X) = names(readPVAR(pfile))
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
