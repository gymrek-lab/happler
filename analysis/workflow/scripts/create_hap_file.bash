#!/usr/bin/env bash

exec 2> "${snakemake_log[0]}"

echo -e "#\torderH\tbeta" > "${snakemake_output[hap]}"
echo -e "#\tversion\t0.1.0" >> "${snakemake_output[hap]}"
echo -e "#H\tbeta\t.2f\tEffect size in linear model" >> "${snakemake_output[hap]}"
echo -e "H\t${snakemake_params[chrom]}\t${snakemake_params[start]}\t${snakemake_params[end]}\thap\t${snakemake_params[beta]}" >> "${snakemake_output[hap]}"

for snp in ${snakemake_input[alleles]}; do
    echo -ne "V\thap\t${snakemake_params[start]}\t${snakemake_params[end]}\t" >> "${snakemake_output[hap]}"
    echo "$snp" | sed 's/\(.*\)-/\1\t/' >> "${snakemake_output[hap]}"
done
