#!/usr/bin/env bash

region="${1:-1_108475598-110283922}"
pop="${2:-EUR_WHITE}"
pheno="${3:-ldl_cholesterol_phenocovar}"

mkdir -p happler_results

happler run \
-t 20 \
--show-tree \
--remove-SNPs \
--out-thresh 20 \
--max-signals 3 \
--chunk-size 500 \
--verbosity DEBUG \
--indep-thresh 15 \
--max-iterations 3 \
-o happler_results/"$region.$pop.hap" \
--discard-multiallelic \
pgens/"$region.$pop.qc.pgen" phenos/"$pheno.residuals.pheno" &>happler_results/"$region.$pop".log

haptools index -o happler_results/"$region.$pop".hap.gz happler_results/"$region.$pop.hap" &>>happler_results/"$region.$pop".log
