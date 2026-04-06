#!/usr/bin/env bash

# Touch all files in a dry-run of Snakemake

set -euo pipefail

LOG_FILE="${1:-log}"
WORKDIR="${WORKDIR:-/tmp/snakemake_gcs_touch}"
mkdir -p "$WORKDIR"

DRY_RUN="${DRY_RUN:-0}"  # DRY_RUN=1 to only print actions

tmpdir="$(mktemp -d "$WORKDIR/touch.XXXXXX")"

# Extract gs:// URIs from output: blocks in appearance order
for uri in $(grep 'output:' "$LOG_FILE" | sed 's/^.*output: //' | sed 's/, /\n/g;s/ \(.* storage\)//g'); do

    # Existence check (skip missing; also skip prefix-style matches)
    ls_out=""
    ls_out="$(gcloud storage ls "$uri" 2>/dev/null || true)"

    if [[ -z "$ls_out" ]]; then
        echo "[SKIP missing] $uri" >&2
        continue
    fi
    if ! grep -Fxq -- "$uri" <<<"$ls_out"; then
        echo "[SKIP non-object/prefix] $uri" >&2
        continue
    fi

    localpath="$tmpdir/object"

    if [[ "$DRY_RUN" == "1" ]]; then
        echo "gcloud storage cp '$uri' '$localpath'"
        echo "gcloud storage cp '$localpath' '$uri'"
    else
        gcloud storage cp "$uri" "$localpath"
        gcloud storage cp "$localpath" "$uri"
    fi

done

rm -rf "$tmpdir"

echo "Done."
