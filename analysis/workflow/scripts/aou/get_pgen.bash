#!/usr/bin/env bash

set -euo pipefail

region="$1"
out_prefix="$2"
threads="${3:-1}"
TEMP_DIR="${4:-}" # default: tmp dir

export GCS_REQUESTER_PAYS_PROJECT="${GOOGLE_PROJECT}"
export GCS_OAUTH_TOKEN="$(gcloud auth application-default print-access-token)"

out_dir="$(dirname "$out_prefix")"
chrom="$(echo "$region" | cut -f1 -d: | sed 's/^chr//')"
pos="$(echo "$region" | cut -f2 -d: | cut -f1 -d-)"
end="$(echo "$region" | cut -f2 -d: | cut -f2 -d-)"
region=chr"${chrom}:${pos}-${end}"
VCF_DIR="${WORKSPACE_BUCKET}/beagle_hg38/chr${chrom}/chr${chrom}."'BATCH*_output.vcf.gz'
CDR_DIR="gs://fc-aou-datasets-controlled/v7"
if [[ -n "$TEMP_DIR" ]]; then
    mkdir -p "$TEMP_DIR"
else
    TEMP_DIR="$(mktemp -d)"
    trap 'rm -rf "$TEMP_DIR"' EXIT
fi
mkdir -p "$TEMP_DIR/$region"

batches="$(gcloud storage ls "$VCF_DIR" | grep -oP '(?<=BATCH)\d+' | sort -n)"
cd "$TEMP_DIR"
echo "workdir: $TEMP_DIR"
echo "$batches" | xargs -P "$threads" -I{} bash -c \
  'batch="$1"; bcftools view --threads 1 -O z -o "$2/$batch.vcf.gz" -r "$3" "$(echo "$4" | sed "s/\*/$batch/")"' _ {} "$region" "$region" "$VCF_DIR"

# assert that all batches were downloaded
for batch in $batches; do
    file_path="$region/$batch".vcf.gz
    if [ ! -f "$file_path" ]; then
        echo "Error: Required file not found at $file_path" >&2
        exit 1
    fi
done

cd -
"$(dirname "$0")"/merge_batched_vcfs.bash "$out_prefix".vcf.gz "$TEMP_DIR/$region"
plink2 --threads "$threads" --out "$out_prefix" --nonfounders --vcf "$out_prefix".vcf.gz --geno 0 --make-pgen --allow-extra-chr --max-alleles 2 --chr "$chrom" --from-bp "$pos" --to-bp "$end"
