#!/usr/bin/env bash

# Efficiently merges VCFs which share variants but have non-overlapping sample sets
# It's important to verify that variants are _exactly_ the same across all VCFs:
# for f in *.vcf.gz; do zcat "$f" | grep -v '^#' | cut -f1-7 | md5sum; done | sort | uniq -c

set -euo pipefail

OUTPUT="$1" # .vcf.gz (can be local path or gs://bucket/path/merged.vcf.gz)
INPUT_DIR="$2" # a dir with a bunch of .vcf.gz files (can be local dir or gs://bucket/dir)

# Helper function to cat files (local or GCS)
cat_file() {
    local f="$1"
    if [[ "$f" == gs://* ]]; then
        gsutil cat "$f"
    else
        cat "$f"
    fi
}

# 1. Collect files using a "version sort" so part2 comes before part10
if [[ "$INPUT_DIR" == gs://* ]]; then
    echo "Detected GCS input. Listing files from bucket..."
    # gsutil ls doesn't support version sort (-v), so we use 'sort -V'
    files=($(gsutil ls "$INPUT_DIR/*.vcf.gz" | sort -V))
else
    echo "Detected local input."
    files=($(ls -v "$INPUT_DIR"/*.vcf.gz))
fi

num_files=${#files[@]}
echo "Found $num_files files. First file is: ${files[0]}"

# 2. Extract the header and first few columns from the first file and then get the GT values from the other files
# Use a subshell (...) to group header output and body output into one stream for bgzip
process_stream() {
    # we turn off pipefail briefly so zcat doesn't kill the script when awk exits early
    set +o pipefail

    # --- A. HEADER PROCESSING ---
    # Logic: 
    # 1. Skip ##INFO lines
    # 2. Skip ##FORMAT lines unless they define ID=GT
    # 3. Print everything else (##fileformat, #CHROM, etc.)
    cat_file "${files[0]}" | zcat | awk '
      /^##INFO/ { next } 
      /^##FORMAT/ { if ($0 ~ /ID=GT/) print; next } 
      /^#CHROM/ { exit } 
      { print }
    '

    set -o pipefail

    # --- B. BODY PROCESSING & MERGING ---
    cmd="paste"

    for ((i=0; i<num_files; i++)); do
        f="${files[$i]}"

        # Build the command string.
        # Note: We must use the specific 'gsutil cat' or 'cat' command inside the process substitution.
        if [[ "$f" == gs://* ]]; then
            # GCS Input
            CAT_CMD="gsutil cat \"$f\""
        else
            # Local Input
            CAT_CMD="cat \"$f\""
        fi

        if [ $i -eq 0 ]; then
            # FILE 1: Columns 1-9 + Samples
            # - $8="."        -> Zap INFO column
            # - sub(/:.*/...) -> Strip everything after ":" in FORMAT (col 9) and Samples (col 10+) to keep only GT
            # - NR==1         : If it's the #CHROM line, print it and move on
            cmd+=" <($CAT_CMD | zcat | grep -v '^##' | awk 'BEGIN{OFS=\"\t\"} NR==1{print; next} {\$8=\".\"; for(i=9;i<=NF;i++) sub(/:.*/, \"\", \$i); print}')"
        else
            # FILES 2-N: Samples Only
            # - cut -f10-     -> Grab sample columns
            # - sed           -> Delete everything after the first colon ":..." until the next tab or end of line to keep only GT
            cmd+=" <($CAT_CMD | zcat | grep -v '^##' | cut -f10- | sed '2,\$s/:[^\\t]*//g')"
        fi
    done

    # --- C. Execute the constructed command ---
    eval "$cmd"
}

# 3. Execute Pipeline
# If OUTPUT is gs://, pipe bgzip output directly to gsutil cp
if [[ "$OUTPUT" == gs://* ]]; then
    echo "Streaming merge directly to GCS: $OUTPUT"
    process_stream | bgzip | gsutil cp - "$OUTPUT"
else
    echo "Writing merge to local file: $OUTPUT"
    process_stream | bgzip > "$OUTPUT"
fi

echo "Done! Output written to $OUTPUT"
