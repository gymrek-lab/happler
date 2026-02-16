#!/usr/bin/env bash

# Efficiently merges VCFs which share variants but have non-overlapping sample sets
# It's important to verify that variants are _exactly_ the same across all VCFs:
# for f in *.vcf.gz; do zcat "$f" | grep -v '^#' | cut -f1-7 | md5sum; done | sort | uniq -c

# This version of the script can resume from where it left. It still needs to be tested.

set -euo pipefail

OUTPUT="$1" # .vcf.gz (can be local path or gs://bucket/path/merged.vcf.gz)
INPUT_DIR="$2" # a dir with a bunch of .vcf.gz files (can be local dir or gs://bucket/dir)
FILTER_IDS="${3:-EnsTR}" # after merging, remove any variants with this pattern in their ID. Set this to "" to disable filtering

# Define number of threads for bgzip (e.g., use all available cores minus 1)
# nproc is a standard linux command to get core count
THREADS=$(nproc)

# Helper function to cat files (local or GCS)
cat_file() {
    local f="$1"
    if [[ "$f" == gs://* ]]; then
        gcloud storage cat "$f"
    else
        cat "$f"
    fi
}

# Helper function to check if file exists (local or GCS)
output_exists() {
    if [[ "$OUTPUT" == gs://* ]]; then
        gcloud storage ls "$OUTPUT" &>/dev/null
    else
        [[ -f "$OUTPUT" ]]
    fi
}

# Check if we're resuming from a partial file
RESUME_FROM=0
if output_exists; then
    echo "Found existing output file. Checking for resumable progress..."
    
    total_lines=$(cat_file "$OUTPUT" | zcat | grep -v '^##' | wc -l || echo 0)

    if [[ "$total_lines" -eq 0 ]]; then
        echo "Partial file appears empty. Delete it first."
        exit 1
    fi

    RESUME_FROM=$((total_lines - 1))
    echo "Resuming from line $total_lines..."
    
    # Rename existing file to .tmp
    if [[ "$OUTPUT" == gs://* ]]; then
        gcloud storage mv "$OUTPUT" "${OUTPUT}.tmp"
    else
        mv "$OUTPUT" "${OUTPUT}.tmp"
    fi
fi

if [[ -n "$FILTER_IDS" ]]; then
    filter_ids() {
        awk -F $'\t' '$3 !~ /'"$FILTER_IDS"'/'
    }
else
    filter_ids() {
        cat
    }
fi

# Skip first RESUME_FROM lines when resuming
skip_lines() {
    if [[ $RESUME_FROM -gt 0 ]]; then
        tail -n +$RESUME_FROM
    else
        cat
    fi
}

# 1. Collect files using a "version sort" so part2 comes before part10
if [[ "$INPUT_DIR" == gs://* ]]; then
    echo "Detected GCS input. Listing files from bucket..."
    # gcloud storage ls doesn't support version sort (-v), so we use 'sort -V'
    files=($(gcloud storage ls "$INPUT_DIR/*.vcf.gz" | sort -V))
else
    echo "Detected local input."
    files=($(ls -v "$INPUT_DIR"/*.vcf.gz))
fi

num_files=${#files[@]}
echo "Found $num_files files. First file is: ${files[0]}"

# 2. Extract the header and first few columns from the first file and then get the GT values from the other files
# Use a subshell (...) to group header output and body output into one stream for bgzip
process_stream() {
    # If resuming, first output the partial file (minus last corrupted line)
    if [[ $RESUME_FROM -gt 0 ]]; then
        echo "Combining first $RESUME_FROM complete lines with resumed output..." >&2
        cat_file "${OUTPUT}.tmp" | zcat | head -n -1
    else
        # Only output header if we're starting from the beginning
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
    fi

    # --- B. BODY PROCESSING & MERGING ---
    cmd="paste"

    for ((i=0; i<num_files; i++)); do
        f="${files[$i]}"

        # Build the command string.
        # Note: We must use the specific 'gcloud storage cat' or 'cat' command inside the process substitution.
        if [[ "$f" == gs://* ]]; then
            # GCS Input
            CAT_CMD="gcloud storage cat \"$f\""
        else
            # Local Input
            CAT_CMD="cat \"$f\""
        fi

        if [ $i -eq 0 ]; then
            # FILE 1: Columns 1-9 + Samples
            # - $8="."        -> Zap INFO column
            # - sub(/:.*/...) -> Strip everything after ":" in FORMAT (col 9) and Samples (col 10+) to keep only GT
            # - NR==1         : If it's the #CHROM line, print it and move on
            cmd+=" <($CAT_CMD | zcat | grep -v '^##' | filter_ids | awk 'BEGIN{OFS=\"\t\"} NR==1{print; next} {\$8=\".\"; for(i=9;i<=NF;i++) sub(/:.*/, \"\", \$i); print}' | skip_lines)"
        else
            # FILES 2-N: Samples Only
            # - cut -f10-     -> Grab sample columns
            # - sed           -> Delete everything after the first colon ":..." until the next tab or end of line to keep only GT
            cmd+=" <($CAT_CMD | zcat | grep -v '^##' | filter_ids | cut -f10- | sed '2,\$s/:[^\\t]*//g' | skip_lines)"
        fi
    done

    # --- C. Execute the constructed command ---
    eval "$cmd"
}

# 3. Execute Pipeline
# If OUTPUT is gs://, pipe bgzip output directly to gcloud storage cp
if [[ "$OUTPUT" == gs://* ]]; then
    echo "Streaming merge directly to GCS: $OUTPUT"
    process_stream | bgzip -@ "$THREADS" | gcloud storage cp - "$OUTPUT"
    [[ $RESUME_FROM -gt 0 ]] && gcloud storage rm "${OUTPUT}.tmp"
else
    echo "Writing merge to local file: $OUTPUT"
    process_stream | bgzip -@ "$THREADS" > "$OUTPUT"
    [[ $RESUME_FROM -gt 0 ]] && rm -f "${OUTPUT}.tmp"
fi

echo "Done! Output written to $OUTPUT"
