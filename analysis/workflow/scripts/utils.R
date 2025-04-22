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
    start_end = as.integer(strsplit(strsplit(region, split=":")[[1]][2], split="-")[[1]])
    matching = pvar$CHROM == chrom & pvar$POS > start_end[1] & pvar$POS < start_end[2]
    pvar = pvar[matching,]
  }
  pos = as.vector(pvar$POS)
  names(pos) = pvar$ID
  pos
}

readPVAR_region = function(pfile, region) {
  chrom = strsplit(region, split=":")[[1]][1]
  start_end = as.integer(strsplit(strsplit(region, split=":")[[1]][2], split="-")[[1]])
  pvar = readPVAR_base(pfile)
  which(pvar$CHROM == chrom & pvar$POS > start_end[1] & pvar$POS < start_end[2])
}

# If provided, the samples argument must be the indices of the desired samples. You can
# retrieve this from the second column of the output of readPSAM()
readPGEN = function(pfile, region=NULL, samples=NULL) {
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
  if (is.null(samples)) {
    X
  } else {
    X[samples,]
  }
}

readPSAM = function(pfile, samples=NULL) {
  if (endsWith(pfile, ".pgen")) {
    pfile = substr(pfile, 1, nchar(pfile)-5)
  }
  psam = data.table::fread(paste0(pfile, '.psam'), sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE)
  if (is.null(samples)) {
    idxs = seq(1, nrow(psam))
  } else {
    idxs = match(intersect(psam[,1], samples), psam[,1])
  }
  data.frame(ids=psam[,1][idxs], idxs)
}

# If a samples subset is specified, the phenotypes will also be reordered to match
readPheno = function(pheno, samples=NULL) {
  pheno = data.table::fread(pheno, sep="\t", header=T, stringsAsFactors = FALSE, check.names=TRUE, data.table=FALSE)
  if (!is.null(samples)) {
    idxs = match(intersect(pheno[,1], samples), pheno[,1])
    pheno = pheno[idxs,]
    # now, reorder so the ordering matches
    idxs = match(samples, pheno[,1])
    pheno = pheno[idxs,]
  }
  pheno
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
  # 
  # Another alternative is to place browser() somewhere in the script, and then run it
  # with the same command line args (but replacing the script name with "R --args").
  # Then, source the script using source() inside the R console

  save(list = ls(all.names = TRUE, envir = env), file = filePath, envir = env)
}

readHap = function(hapfile) {
  # Return a data.frame with the start, end, and ID columns of the haplotype lines
  # in the .hap file
  # Parameters:
  #   hapfile: The path to the .hap file
  data <- tryCatch(read.table(hapfile, header = FALSE, fill = TRUE, stringsAsFactors = FALSE), error=function(e) data.frame())
  if (nrow(data) != 0) {
    # Filter rows that start with "H" and select columns 3 to 5
    h_data <- data[grep("^H", data$V1), c(3, 4, 5)]
    # Rename the columns for clarity (optional)
    colnames(h_data) <- c("start", "end", "id")
    h_data
  } else {
    data
  }
}
