#!/usr/bin/env bash

# Touch all files in a dry-run of Snakemake

set -euo pipefail

LOG_FILE="${1:log}"
WORKDIR="${WORKDIR:-/tmp/snakemake_gcs_touch}"
mkdir -p "$WORKDIR"

DRY_RUN="${DRY_RUN:-0}"  # DRY_RUN=1 to only print actions

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing required command: $1" >&2; exit 1; }; }
need_cmd gcloud
need_cmd awk
need_cmd sed
need_cmd mktemp
need_cmd grep

# Extract gs:// URIs from output: blocks, in appearance order, and de-dupe while preserving order.
mapfile -t GCS_URIS < <(
  awk '
    BEGIN {inout=0}
    /^    output:/ {inout=1; sub(/^    output:[[:space:]]*/, ""); print; next}
    inout==1 {
      if ($0 ~ /^    /) { sub(/^    /,""); print; next }
      inout=0
    }
  ' "$LOG_FILE" \
  | tr ',' '\n' \
  | sed -E 's/\(send to storage\)//g; s/\(retrieve from storage\)//g; s/^[[:space:]]+//; s/[[:space:]]+$//' \
  | awk '
      /^gs:\/\// {
        if (!seen[$0]++) print $0
      }
    '
)

if [[ ${#GCS_URIS[@]} -eq 0 ]]; then
  echo "No gs:// URIs found in output: blocks in $LOG_FILE" >&2
  exit 2
fi

echo "Found ${#GCS_URIS[@]} unique output objects (order preserved)."
echo "Workdir: $WORKDIR"
echo "Dry-run: $DRY_RUN"
echo

touch_one() {
  local uri="$1"

  # Existence check (skip missing; also skip prefix-style matches)
  local ls_out=""
  ls_out="$(gcloud storage ls "$uri" 2>/dev/null || true)"

  if [[ -z "$ls_out" ]]; then
    echo "[SKIP missing] $uri"
    return 0
  fi
  if ! grep -Fxq -- "$uri" <<<"$ls_out"; then
    echo "[SKIP non-object/prefix] $uri"
    return 0
  fi

  local tmpdir localpath
  tmpdir="$(mktemp -d "$WORKDIR/touch.XXXXXX")"
  localpath="$tmpdir/object"

  if [[ "$DRY_RUN" == "1" ]]; then
    echo "[DRY] gcloud storage cp '$uri' '$localpath'"
    echo "[DRY] gcloud storage cp '$localpath' '$uri'"
  else
    gcloud storage cp "$uri" "$localpath"
    gcloud storage cp "$localpath" "$uri"
  fi

  rm -rf "$tmpdir"
  echo "[OK] touched $uri"
}

# Sequential, in-order processing:
for uri in "${GCS_URIS[@]}"; do
  touch_one "$uri"
done

echo
echo "Done."
