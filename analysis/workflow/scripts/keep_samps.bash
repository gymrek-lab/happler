#!/usr/bin/env bash

# arg1: VCF file from which to obtain SNP genotypes
# arg2: VCF file from which to obtain STR genotypes
# arg3: a list of samples to exclude
# ex: workflow/scripts/keep_samps.bash out/1_98001984-99001984/phased.bcf strs/chr1.vcf.gz data/ukb_withdrawals.tsv | head




snp_vcf="$1"
str_vcf="$2"
samps_file="$3"


comm -13 \
    <(sort "$samps_file") \
    <(bcftools query -l "$str_vcf" | cut -d'_' -f1 | grep -vE '^\-' | sort) | \
comm -12 - \
    <(bcftools query -l "$snp_vcf" | sort)
