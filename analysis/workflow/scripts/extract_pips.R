#!/usr/bin/env Rscript

# This R script extracts PIPs from a SuSiE RDS file produced by run_SuSiE.R

# param1: An RDS file output by run_SuSiE.R
# params2: A TSV file to which to write the PIPs

args = commandArgs(trailingOnly = TRUE)
susie_rds = args[1]
out = args[2]

susie_results = readRDS(susie_rds)
write.table(susieR::susie_get_pip(susie_results$fitted), file=out, row.names=T, sep="\t", quote=F, col.names=F)
