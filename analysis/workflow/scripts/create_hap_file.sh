#!/usr/bin/env bash

exec 2> "${snakemake_log}"

echo -e "#\torderH\tbeta" > "${snakemake_output[hap]}"
echo -e "#\tversion\t0.1.0" >> "${snakemake_output[hap]}"
echo -e "#H\tbeta\t.2f\tEffect size in linear model" >> "${snakemake_output[hap]}"

starts=()
ends=()
for snp in ${snakemake_params[4]}; do
    snp_id="$(echo "$snp" | sed 's/\(.*\):/\1\t/')"
    start_pos="$(grep -m1 "$(echo "$snp_id" | cut -f1)" "${snakemake_input[gts]}" | cut -f2)"
    end_pos=$((start_pos+1))
    echo -e "V\thap\t$start_pos\t$end_pos\t$snp_id" >> "${snakemake_output[hap]}"
    starts+=("$start_pos")
    ends+=("$end_pos")
done

start_pos="$(echo "${starts[@]}" | tr ' ' '\n' | sort -n | head -n1)"
end_pos="$(echo "${ends[@]}" | tr ' ' '\n' | sort -rn | head -n1)"
echo -e "H\t${snakemake_params[1]}\t$start_pos\t$end_pos\thap\t${snakemake_params[3]}" >> "${snakemake_output[hap]}"

LC_ALL=C sort -k1,4 -o "${snakemake_output[hap]}" "${snakemake_output[hap]}"
