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
    # A. Print Header (removing ##INFO lines)
    zcat "${files[0]}" | grep "^#" | grep -v "^##INFO"

    # B. Construct the dynamic paste command
    cmd="paste"

    for ((i=0; i<num_files; i++)); do
        f="${files[$i]}"
        
        if [ $i -eq 0 ]; then
            # FIRST FILE: Remove header, set INFO (col 8) to ".", keep all columns
            # Note: We escape $8 as \$8 so it isn't evaluated by the shell now
            cmd+=" <(zcat \"$f\" | grep -v '^#' | awk 'BEGIN{OFS=\"\t\"} {\$8=\".\"; print}')"
        else
            # OTHER FILES: Remove header, cut columns 10-End
            cmd+=" <(zcat \"$f\" | grep -v '^#' | cut -f10-)"
        fi
    done

    # C. Execute the constructed command
    eval "$cmd"

) | bgzip > "$OUTPUT"

echo "Done! Output written to $OUTPUT"
