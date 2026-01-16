#!/usr/bin/env bash

set -euo pipefail

OUTPUT="$1" # .vcf.gz
in_vcfs_dir="$2" # a dir with a bunch of .vcf.gz files

# 1. Collect files using a "Version Sort" so part2 comes before part10
files=($(ls -v "$in_vcfs_dir"/*.vcf.gz))
num_files=${#files[@]}

echo "Found $num_files files. First file is: ${files[0]}"

# 2. Extract and clean the Header from the first file
# Use a subshell (...) to group header output and body output into one stream for bgzip
(
    # --- A. HEADER PROCESSING ---
    # Logic: 
    # 1. Skip ##INFO lines
    # 2. Skip ##FORMAT lines unless they define ID=GT
    # 3. Print everything else (##fileformat, #CHROM, etc.)
    zcat "${files[0]}" | grep "^#" | \
    awk '
      /^##INFO/ { next } 
      /^##FORMAT/ { if ($0 ~ /ID=GT/) print; next } 
      { print }
    '

    # --- B. BODY PROCESSING & MERGING ---
    cmd="paste"

    for ((i=0; i<num_files; i++)); do
        f="${files[$i]}"
        
        if [ $i -eq 0 ]; then
            # FILE 1: Columns 1-9 + Samples
            # - $8="."        -> Zap INFO column
            # - sub(/:.*/...) -> Strip everything after ":" in FORMAT (col 9) and Samples (col 10+) to keep only GT
            cmd+=" <(zcat \"$f\" | grep -v '^#' | awk -F $'\t' 'BEGIN{OFS=\"\t\"} {\$8=\".\"; for(i=9;i<=NF;i++) sub(/:.*/, \"\", \$i); print}')"
        else
            # FILES 2-N: Samples Only
            # - cut -f10-     -> Grab sample columns
            # - sub(/:.*/...) -> Strip everything after ":" in all columns to keep only GT
            cmd+=" <(zcat \"$f\" | grep -v '^#' | cut -f10- | awk 'BEGIN{OFS=\"\t\"} {for(i=1;i<=NF;i++) sub(/:.*/, \"\", \$i); print}')"
        fi
    done

    # C. Execute the constructed command
    eval "$cmd"

) | bgzip > "$OUTPUT"

echo "Done! Output written to $OUTPUT"
