#!/usr/bin/env bash

set -euo pipefail

OUTPUT="$1" # .vcf.gz
in_vcfs_dir="$2" # a dir with a bunch of .vcf.gz files

# 1. Collect files using a "version sort" so part2 comes before part10
files=($(ls -v "$in_vcfs_dir"/*.vcf.gz))
num_files=${#files[@]}

echo "Found $num_files files. First file is: ${files[0]}"

# 2. Extract and clean the Header from the first file
# Use a subshell (...) to group header output and body output into one stream for bgzip
(
    # We turn off pipefail briefly so zcat doesn't kill the script when awk exits early
    set +o pipefail

    # --- A. HEADER PROCESSING ---
    # Logic: 
    # 1. Skip ##INFO lines
    # 2. Skip ##FORMAT lines unless they define ID=GT
    # 3. Print everything else (##fileformat, #CHROM, etc.)
    zcat "${files[0]}" | awk '/^##INFO/ { next } /^##FORMAT/ { if ($0 ~ /ID=GT/) print; next } /^#CHROM/ { exit } { print }'

    set -o pipefail # Turn safety back on

    # --- B. BODY PROCESSING & MERGING ---
    cmd="paste"

    for ((i=0; i<num_files; i++)); do
        f="${files[$i]}"
        
        if [ $i -eq 0 ]; then
            # FILE 1: Columns 1-9 + Samples
            # - $8="."        -> Zap INFO column
            # - sub(/:.*/...) -> Strip everything after ":" in FORMAT (col 9) and Samples (col 10+) to keep only GT
            # - NR==1         : If it's the #CHROM line, print it and move on
            cmd+=" <(zcat \"$f\" | grep -v '^##' | awk 'BEGIN{OFS=\"\t\"} NR==1{print; next} {\$8=\".\"; for(i=9;i<=NF;i++) sub(/:.*/, \"\", \$i); print}')"
        else
            # FILES 2-N: Samples Only
            # - cut -f10-     -> Grab sample columns
            # - sed           -> Delete everything after the first colon ":..." until the next tab or end of line to keep only GT
            cmd+=" <(zcat \"$f\" | grep -v '^##' | cut -f10- | sed '2,\$s/:[^\\t]*//g')"
        fi
    done

    # C. Execute the constructed command
    eval "$cmd"

) | bgzip > "$OUTPUT"

echo "Done! Output written to $OUTPUT"
