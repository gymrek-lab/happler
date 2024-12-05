#!/usr/bin/env bash

exec 2> "${snakemake_log}"

echo -e "#\tversion\t0.2.0" >> "${snakemake_output[hap]}"
echo -e "#R\tbeta\t.2f\tEffect size in linear model" >> "${snakemake_output[hap]}"

if [ "${snakemake_params[repeat]}" -eq "0" ]; then
    starts=()
    ends=()
    for snp in ${snakemake_params[4]}; do
        snp_id="$(echo "$snp" | sed 's/\(.*\):/\1\t/')"
        start_pos="$(grep -Ev '^#' -m1 "$(echo "$snp_id" | cut -f1)" "${snakemake_input[gts]}" | cut -f2)"
        end_pos=$((start_pos+1))
        echo -e "V\thap\t$start_pos\t$end_pos\t$snp_id" >> "${snakemake_output[hap]}"
        starts+=("$start_pos")
        ends+=("$end_pos")
    done

    start_pos="$(echo "${starts[@]}" | tr ' ' '\n' | sort -n | head -n1)"
    end_pos="$(echo "${ends[@]}" | tr ' ' '\n' | sort -rn | head -n1)"
    echo -e "H\t${snakemake_params[1]}\t$start_pos\t$end_pos\thap\t${snakemake_params[3]}" >> "${snakemake_output[hap]}"
else
    for str in ${snakemake_params[5]}; do
        str_id=$'\t'"$str"$'\t'
        start_pos="$(grep --line-buffered -Ev '^#' "${snakemake_input[gts]}" | grep -Pm1 "${start_pos}${str_id}" | cut -f2)"
        end_pos="$(grep --line-buffered -Ev '^#' "${snakemake_input[gts]}" | grep -Pm1 "${start_pos}${str_id}" | grep -Po '(?<=;END=)\d+(?=;PERIOD)')"
        echo -e "R\t${snakemake_params[1]}\t$start_pos\t$end_pos\t$str\t${snakemake_params[3]}" >> "${snakemake_output[hap]}"
    done
fi

LC_ALL=C sort -k1,4 -o "${snakemake_output[hap]}" "${snakemake_output[hap]}"
