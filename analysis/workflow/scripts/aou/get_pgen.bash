#!/usr/bin/env bash

set -euo pipefail

region="$1"
out_prefix="$2"

export GCS_REQUESTER_PAYS_PROJECT="${GOOGLE_PROJECT}"
export GCS_OAUTH_TOKEN="$(gcloud auth application-default print-access-token)"

out_dir="$(dirname "$out_prefix")"
chrom="$(echo "$region" | cut -f1 -d: | sed 's/chr//')"
pos="$(echo "$region" | cut -f2 -d: | cut -f1 -d-)"
end="$(echo "$region" | cut -f2 -d: | cut -f2 -d-)"
VCF_DIR="${WORKSPACE_BUCKET}/beagle_hg38/chr${chrom}/chr${chrom}."'BATCH*_output.vcf.gz'
CDR_DIR="gs://fc-aou-datasets-controlled/v7"

batches="$(gsutil ls "$VCF_DIR" | grep -oP '(?<=BATCH)\d+' | sort -n)"
cd "$out_dir/batches"
for batch in $batches; do
    bcftools view -O b -o "$batch".bcf -r "$region" "$(echo "$VCF_DIR" | sed 's/*/'"$batch"'/')"
    # plink2 --out "$batch" --nonfounders --bcf "$batch".bcf --geno 0 --make-pgen --allow-extra-chr --max-alleles 2 --chr "$chrom" --from-bp "$pos" --to-bp "$end"
done

# assert that all batches were downloaded
for batch in $batches; do
    file_path="$batch".bcf
    if [ ! -f "$file_path" ]; then
        echo "Error: Required file not found at $file_path" >&2
        exit 1
    fi
fi

cd -
bcftools merge --no-index -O b -o "$out_prefix".bcf -l <(ls "$out_dir/batches"/*.bcf)
plink2 --out "$out_prefix" --nonfounders --bcf "$out_prefix".bcf --geno 0 --make-pgen --allow-extra-chr --max-alleles 2 --chr "$chrom" --from-bp "$pos" --to-bp "$end"
gsutil cp "$out_prefix".p{gen,var,sam} ${WORKSPACE_BUCKET}/aryarm/pgens/ALL_SAMPLES/
