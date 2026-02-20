#!/usr/bin/env bash

# Efficiently merges VCFs which share variants but have non-overlapping sample sets
# It's important to verify that variants are _exactly_ the same across all VCFs:
# for f in *.vcf.gz; do zcat "$f" | grep -v '^#' | cut -f1-7 | md5sum; done | sort | uniq -c

# Note that the script can resume where it left off if it gets interrupted.

# To test and benchmark this script, you can download the first few variants of all batches in chr22 and try to merge them with this script vs bcftools.
# Then, convert them to PGEN and compare them with plink2 --pgen-diff to make sure they are the same.
# TODO: the code for that

set -euo pipefail

OUTPUT="$1" # .vcf.gz (can be local path or gs://bucket/path/merged.vcf.gz)
INPUT_DIR="$2" # a dir with a bunch of .vcf.gz files (can be local dir or gs://bucket/dir)
FILTER_IDS="${3:-EnsTR}" # after merging, remove any variants with this pattern in their ID. Set this to "" to disable filtering

# Define number of threads for bgzip (e.g., use all available cores minus 1)
# nproc is a standard linux command to get core count
THREADS=$(nproc)

# Helper function to cat files (local or GCS)
read_file() {
    local f="$1"
    if [[ "$f" == gs://* ]]; then
        gcloud storage cat "$f"
    else
        cat "$f"
    fi
}

write_file() {
    local f="$1"
    if [[ "$f" == gs://* ]]; then
        gcloud storage cp - "$f"
    else
        cat > "$f"
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
RESUME_REGION=""
if output_exists; then
    echo "Found existing output file. Assuming previously interrupted. Resuming..."
    export GCS_OAUTH_TOKEN="$(gcloud auth application-default print-access-token)"

    # Extract the genomic position of the last complete variant line
    # Also, copy all but the last line (which is likely truncated) to a temporary file
    last_variant="$({ { read_file "${OUTPUT}" | zcat; } 2>/dev/null || true; } | tee >(head -n -1 | bgzip -@ "$THREADS" | write_file "${OUTPUT}.tmp") | tail -n 1)"

    if [[ -z "$last_variant" && "$last_variant" != "#"* ]]; then
        echo "Found existing output but file appears empty. Delete it first."
        rm -f "${OUTPUT}.tmp"
        exit 1
    fi

    RESUME_REGION="$(echo "$last_variant" | cut -f1):$(echo "$last_variant" | cut -f2)"
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

# 1. Collect files using a "version sort" so part2 comes before part10
if [[ "$INPUT_DIR" == gs://* ]]; then
    echo "Detected GCS input"
    # gcloud storage ls doesn't support version sort (-v), so we use 'sort -V'
    files=($(gcloud storage ls "$INPUT_DIR/*.vcf.gz" | sort -V))
else
    echo "Detected local input"
    files=($(ls -v "$INPUT_DIR"/*.vcf.gz))
fi

num_files=${#files[@]}
echo "Found $num_files files. First file is: ${files[0]}"

# 2. Extract the header and first few columns from the first file and then get the GT values from the other files
# Use a subshell (...) to group header output and body output into one stream for bgzip
process_stream() {
    # Only output header if we're not resuming an interrupted file
    if [[ -z "$RESUME_REGION" ]]; then
        # we turn off pipefail briefly so zcat doesn't kill the script when awk exits early
        set +o pipefail

        # --- A. HEADER PROCESSING ---
        # Logic: 
        # 1. Skip ##INFO lines
        # 2. Skip ##FORMAT lines unless they define ID=GT
        # 3. Print everything else (##fileformat, #CHROM, etc.)
        read_file "${files[0]}" | zcat | awk '
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

        if [ $i -eq 0 ]; then
            # FILE 1: Columns 1-9 + Samples
            # - $8="."        -> Zap INFO column
            # - sub(/:.*/...) -> Strip everything after ":" in FORMAT (col 9) and Samples (col 10+) to keep only GT
            # - NR==1         : If it's the #CHROM line, print it and move on (only when not resuming)
            if [[ -z "$RESUME_REGION" ]]; then
                cmd+=" <(read_file \"$f\" | zcat | grep -v '^##' | filter_ids | awk 'BEGIN{OFS=\"\t\"} NR==1{print; next} {\$8=\".\"; for(i=9;i<=NF;i++) sub(/:.*/, \"\", \$i); print}')"
            else
                # If resuming an interrupted file, use tabix to skip lines. Also, no need to process the header anymore.
                cmd+=" <(tabix \"$f\" \"${RESUME_REGION}-\" | filter_ids | awk 'BEGIN{OFS=\"\t\"} {\$8=\".\"; for(i=9;i<=NF;i++) sub(/:.*/, \"\", \$i); print}')"
            fi
        else
            # FILES 2-N: Samples Only
            # - cut -f10-     -> Grab sample columns
            # - sed           -> Delete everything after the first colon ":..." until the next tab or end of line to keep only GT
            # - sed 2,        -> Skip processing of the first line (#CHROM)
            if [[ -z "$RESUME_REGION" ]]; then
                cmd+=" <(read_file \"$f\" | zcat | grep -v '^##' | filter_ids | cut -f10- | sed '2,\$s/:[^\\t]*//g')"
            else
                # If resuming an interrupted file, use tabix to skip lines. Also, no need to process the header anymore.
                cmd+=" <(tabix \"$f\" \"${RESUME_REGION}-\" | filter_ids | cut -f10- | sed 's/:[^\\t]*//g')"
            fi
        fi
    done

    # --- C. Execute the constructed command ---
    eval "$cmd"
}

# 3. Execute Pipeline
# If OUTPUT is gs://, pipe bgzip output directly to gcloud storage cp
if [[ "$OUTPUT" == gs://* ]]; then
    echo "Streaming merge directly to GCS: $OUTPUT"
    if [[ -n "$RESUME_REGION" ]]; then
        # We can't append to an existing file in GCS
        read_file "${OUTPUT}.tmp"
        gcloud storage rm "${OUTPUT}.tmp"
        process_stream | bgzip -@ "$THREADS"
    else
        process_stream | bgzip -@ "$THREADS"
    fi | write_file "$OUTPUT"
else
    echo "Writing merge to local file: $OUTPUT"
    if [[ -n "$RESUME_REGION" ]]; then
        mv "${OUTPUT}.tmp" "${OUTPUT}"
        process_stream | bgzip -@ "$THREADS" >> "${OUTPUT}"
    else
        process_stream | bgzip -@ "$THREADS" | write_file "${OUTPUT}"
    fi
fi

echo "Done! Output written to $OUTPUT"
