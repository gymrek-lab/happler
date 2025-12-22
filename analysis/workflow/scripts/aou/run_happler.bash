#!/usr/bin/env bash

region="${1:-1_108475598-110283922}"
pop="${2:-EUR_WHITE}"
pheno="${3:-ldl_cholesterol_phenocovar}"

mkdir -p happler
cd happler

happler run \
-t 20 \
--show-tree \
--remove-SNPs \
--out-thresh 20 \
--max-signals 3 \
--verbosity DEBUG \
--indep-thresh 15 \
--max-iterations 3 \
-o "$region.$pop.hap" \
--discard-multiallelic \
pgens/"$region.$pop.pgen" phenos/"$pheno.residuals.pheno" &>"$region.$pop".log

haptools index -o "$region.$pop".hap.gz "$region.$pop.hap"
