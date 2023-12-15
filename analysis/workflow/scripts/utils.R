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
