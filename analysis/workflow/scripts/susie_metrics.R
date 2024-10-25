#!/usr/bin/env Rscript

# This R script computes specific metrics from the output of a SuSiE run

# param1: The path to a RDS file containing the output of SuSiE
# param2: The path to a .hap file containing the observed haplotype

# Output will be written in tab-separated format to stdout
# To debug this script, add a browser() call somewhere in it. Then run it with its command line args but replace the script name with "R --args".

# first, we import the args
args = commandArgs(trailingOnly = TRUE)
susie_res = args[1]
obs_hap = args[2]

# load functions to help read various file types
source("workflow/scripts/utils.R")

write("Parsing observed hap file", stderr())
happler_hap = readHap(obs_hap)

write("Parsing SuSiE results", stderr())
# parse the susie data:
# 1) fitted - the list provided by susie as output
# 2) susie_pip - the PIPs output by SuSiE
susie_res = readRDS(susie_res)
fitted = susie_res$fitted
susie_pip = fitted$pip

for (happler_hap_file_idx in 1:nrow(happler_hap)) {

    happler_hap_id = happler_hap[happler_hap_file_idx, "id"]
    happler_hap_idx = which(names(susie_pip) %in% happler_hap_id)

    write(paste("Computing metrics for", happler_hap_id), stderr())
    # The metrics are:
    # 1) What is the PIP of the observed hap?
    # 2) Does the observed hap get the highest PIP?
    # 3) What is the best PIP among the variants?
    # 4) Is the observed hap in a credible set?
    # 5) What is the purity of the credible set?
    # 6) What is the length of the credible set?
    obs_pip = susie_pip[happler_hap_id]
    has_highest_pip = as.integer(sum(obs_pip < susie_pip) == 0)
    # TODO: fix this so that it also excludes other haplotypes besides obs_pip
    best_variant_pip = max(susie_pip[names(susie_pip) != names(obs_pip)])
    credible_set = NULL
    for (cs in names(fitted$sets$cs)) {
        if (happler_hap_idx %in% fitted$sets$cs[cs]) {
            credible_set = cs
            break
        }
    }
    purity = 0
    cs_length = 0
    if (!is.null(credible_set)) {
        purity = fitted$sets$purity[cs, "mean.abs.corr"]
        cs_length = length(fitted$sets$cs[credible_set])
    }
    in_credible_set = as.integer(!is.null(credible_set))

    write(paste("Outputting metrics for", happler_hap_id), stderr())
    write(paste(obs_pip, has_highest_pip, best_variant_pip, in_credible_set, purity, cs_length), stdout())

}
