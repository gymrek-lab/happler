#!/usr/bin/env bash

exec 2> "${snakemake_log}"

echo -e "#\torderH\tbeta" > "${snakemake_output[hap]}"
echo -e "#\tversion\t0.1.0" >> "${snakemake_output[hap]}"
echo -e "#H\tbeta\t.2f\tEffect size in linear model" >> "${snakemake_output[hap]}"
echo -e "H\t${snakemake_params[1]}\t${snakemake_params[2]}\thap\t${snakemake_params[3]}" >> "${snakemake_output[hap]}"

for snp in ${snakemake_params[4]}; do
    # TODO: figure out how to use the proper start and end positions
    echo -ne "V\thap\t${snakemake_params[2]}\t" >> "${snakemake_output[hap]}"
    echo "$snp" | sed 's/\(.*\):/\1\t/' >> "${snakemake_output[hap]}"
done
