#!/usr/bin/env Rscript

# This R script computes specific metrics from the output of a SuSiE run

# param1: The path to a RDS file containing the output of SuSiE
# param2: The path to a .hap file containing the observed haplotype
# param3: The path to a .hap file containing the causal haplotype

# Output will be written in tab-separated format to stdout
# To debug this script, add a browser() call somewhere in it. Then run it with its command line args but replace the script name with "R --args".

# first, we import the args
args = commandArgs(trailingOnly = TRUE)
susie_res = args[1]
obs_hap = args[2]
caus_hap = args[3]

write("Parsing observed hap file", stderr())
for (line in readLines(obs_hap)) {
    # we assume the first line in the hap file that begins with H denotes the haplotype
    if (grepl("^H\\t", line)) break
}
# extract the start, end, and ID of the haplotype
happler_hap = as.integer(strsplit(line, "\t")[[1]][c(3,4)])
happler_hap_id = strsplit(line, "\t")[[1]][c(5)]

write("Parsing causal hap file", stderr())
for (line in readLines(caus_hap)) {
    # we assume the first line in the hap file that begins with H denotes the haplotype
    if (grepl("^H\\t", line)) break
}
# extract the start, end, and ID of the haplotype
causal_hap = as.integer(strsplit(line, "\t")[[1]][c(3,4)])
causal_hap_id = strsplit(line, "\t")[[1]][c(5)]

write("Parsing SuSiE results", stderr())
# parse the susie data:
# 1) fitted - the list provided by susie as output
# 2) susie_pip - the PIPs output by SuSiE
susie_res = readRDS(susie_res)
fitted = susie_res$fitted
susie_pip = fitted$pip
happler_hap_idx = which(names(susie_pip) %in% happler_hap_id)

write("Computing metrics", stderr())
# The metrics are:
# 1) What is the PIP of the observed hap?
# 2) Does the observed hap get the highest PIP?
# 3) What is the best PIP among the variants?
# 4) Is the observed hap in a credible set?
# 5) What is the purity of the credible set?
obs_pip = susie_pip[happler_hap_id]
has_highest_pip = as.integer(sum(obs_pip < susie_pip) == 0)
best_variant_pip = max(susie_pip[names(susie_pip) != names(obs_pip)])
credible_set = NULL
for (cs in names(fitted$sets$cs)) {
    if (happler_hap_idx %in% fitted$sets$cs[cs]) {
        credible_set = cs
        break
    }
}
purity = 0
if (!is.null(credible_set)) {
    purity = fitted$sets$purity[cs, "mean.abs.corr"]
}
in_credible_set = as.integer(!is.null(credible_set))

write("Outputting metrics", stderr())
write(paste(obs_pip, has_highest_pip, best_variant_pip, in_credible_set, purity), stdout())
