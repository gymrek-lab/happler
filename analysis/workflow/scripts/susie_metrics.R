#!/usr/bin/env Rscript

# This R script computes specific metrics from the output of a SuSiE run

# param1: The path to a RDS file containing the output of SuSiE
# param2: The path to a .hap file containing the observed/causal haplotype

# Output will be written in tab-separated format to stdout
# To debug this script, add a browser() call somewhere in it. Then run it with its command line args (but replacing the script name with "R --args") and then source it


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

# first, we import the args
args = commandArgs(trailingOnly = TRUE)
gt = args[1]
phen = args[2]
susie_res = args[3]
curr_hap = args[4]
exclude_causal = args[5]
region = args[6]

# load functions to help read various file types
source(paste0(thisDir, "/utils.R"))

write("Reading sample list and phenotypes", stderr())
samples = readPSAM(gt, samples=readPheno(phen)[,1])
phen = readPheno(phen, samples=samples[,1])
# load psam and ensure sample names are the same
stopifnot(samples[,1] == phen[,1])
y = as.matrix(phen[,2])
# import genotype matrices as proper matrices
write(paste("Trying to read genotype matrix with", nrow(samples), "samples"), stderr())
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

happler_hap_ids = NULL
# TODO: also handle case where it's a snplist file instead
if (curr_hap != "NULL") {
    write("Parsing observed/causal hap file", stderr())
    if (endsWith(curr_hap, ".hap")) {
        happler_hap_ids = readHap(curr_hap)[,"id"]
    } else {
        # if it's not a hap file, we assume it's a snplist file
        happler_hap_ids = read.table(curr_hap, header = FALSE, stringsAsFactors = FALSE)[,1]
    }
}

write("Parsing SuSiE results", stderr())
# parse the susie data:
# 1) fitted - the list provided by susie as output
# 2) susie_pip - the PIPs output by SuSiE
susie_res = readRDS(susie_res)
fitted = susie_res$fitted
susie_pip = susieR::susie_get_pip(fitted)
names(susie_pip) = names(fitted$pip)

# if a haplotype was not provided, we just output metrics for each credible set
if (is.null(happler_hap_ids)) {
    # do not filter susie credible sets by default purity or coverage
    susie_CS_details = susieR::susie_get_cs(fitted, X = X, min_abs_corr = 0)
    susie_CSs = susie_CS_details$cs
    write("No haplotype provided", stderr())
    # The metrics are:
    # 1) False (NA)
    # 2) False (0)
    # 3) False (0)
    # 4) What is the best PIP in the credible set?
    # 5) False (0)
    # 6) How many credible sets are there?
    # 7) What is the purity of this credible set?
    # 8) What is the length of this credible set?
    obs_pip = 0
    coverage = 0
    curr_id = "NA"
    has_highest_pip = 0
    in_credible_set = 0
    num_credible_sets = length(susie_CSs)
    if (length(susie_CSs) == 0) {
        # if there are no credible sets, we warn the user
        write("No credible sets found", stderr())
    } else {
        # iterate through each of the credible sets to find the one with the best variant in it
        for (credible_set in names(susie_CSs)) {
            best_variant_pip = max(susie_pip[susie_CSs[[credible_set]]])
            purity = susie_CS_details$purity[credible_set, "min.abs.corr"]
            cs_length = length(susie_CSs[[credible_set]])
            write(paste(curr_id, obs_pip, has_highest_pip, best_variant_pip, in_credible_set, num_credible_sets, purity, cs_length), stdout())   
        }
    }
    quit()
}

# filter susie credible sets by default purity
susie_CS_details = susieR::susie_get_cs(fitted, X = X)
susie_CSs = susie_CS_details$cs

for (happler_hap_id in happler_hap_ids) {

    happler_hap_idx = which(names(susie_pip) %in% happler_hap_id)

    write(paste("Computing metrics for", happler_hap_id), stderr())
    # The metrics are:
    # 1) What is the ID of the hap?
    # 2) What is the PIP of the hap?
    # 3) Does the observed/causal hap get the highest PIP?
    # 4) What is the next best PIP in the credible set, excluding the hap?
    # 5) Is the observed/causal hap in a credible set? If so, what is it's index?
    # 6) How many credible sets are there?
    # 7) What is the purity of the credible set with the observed/causal hap?
    # 8) What is the length of the credible set with the observed/causal hap?
    obs_pip = susie_pip[happler_hap_id]
    has_highest_pip = as.integer(sum(obs_pip < susie_pip) == 0)
    credible_set = NULL
    num_credible_sets = length(susie_CSs)
    # iterate through each of the credible sets to find the one with the hap in it
    for (cs in names(susie_CSs)) {
        if (happler_hap_idx %in% susie_CSs[[cs]]) {
            credible_set = cs
            break
        }
    }
    purity = 0
    coverage = 0
    cs_length = 0
    in_credible_set = 0
    best_variant_pip = 0
    if (!is.null(credible_set)) {
        purity = susie_CS_details$purity[credible_set, "min.abs.corr"]
        coverage = susie_CS_details$purity[credible_set, "coverage"]
        actual_cs = susie_CSs[[credible_set]]
        cs_length = length(actual_cs)
        # if the hap was in a credible set, what is the index of the credible set?
        in_credible_set = which(names(susie_CSs) == credible_set)
        # what was the next best PIP in the credible set, excluding the hap?
        if (cs_length > 1) {
            # if the credible set has more than one variant, we can find the next best PIP
            other_variant_idxs = actual_cs[!happler_hap_idx==actual_cs]
            best_variant_pip = max(susie_pip[other_variant_idxs])
        }
    }

    write(paste("Outputting metrics for", happler_hap_id), stderr())
    write(paste(happler_hap_id, obs_pip, has_highest_pip, best_variant_pip, in_credible_set, num_credible_sets, purity, cs_length), stdout())

}
