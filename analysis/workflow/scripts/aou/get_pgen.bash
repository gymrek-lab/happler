#!/usr/bin/env bash

region="$1"
pheno="$2"

export GCS_REQUESTER_PAYS_PROJECT="${GOOGLE_PROJECT}"
export GCS_OAUTH_TOKEN="$(gcloud auth application-default print-access-token)"

out_prefix="$(echo "$region" | sed 's/:/_/;s/chr//')"
mkdir -p "$out_prefix"
chrom="$(echo "$region" | cut -f1 -d: | sed 's/chr//')"
pos="$(echo "$region" | cut -f2 -d: | cut -f1 -d-)"
end="$(echo "$region" | cut -f2 -d: | cut -f2 -d-)"
VCF_DIR="${WORKSPACE_BUCKET}/beagle_hg38/chr${chrom}/chr${chrom}."'BATCH*_output.vcf.gz'
CDR_DIR="gs://fc-aou-datasets-controlled/v7"

batches="$(gsutil ls "$VCF_DIR" | grep -oP '(?<=BATCH)\d+' | sort -n)"
cd "$out_prefix"
for batch in $batches; do
    bcftools view -O b -o "$batch".bcf -r "$region" "$(echo "$VCF_DIR" | sed 's/*/'"$batch"'/')"
#    plink2 --out "$batch" --nonfounders --bcf "$batch".bcf --geno 0 --make-pgen --allow-extra-chr --max-alleles 2 --chr "$chrom" --from-bp "$pos" --to-bp "$end"
done

cd ..   
bcftools merge --no-index -O z -o "$out_prefix".vcf.bgz -l <(ls "$out_prefix"/*.bcf)

# now, let's use hail to filter the GT data
../happler/analysis/workflow/scripts/aou/hail_qc_EUR_AFR.py "$out_prefix".vcf.bgz "$pheno"

# note that we skip --maf bc the input is already filtered
plink2 --out "$out_prefix" --nonfounders --vcf "$out_prefix".qc.vcf.bgz --geno 0 --make-pgen --allow-extra-chr --max-alleles 2 --chr "$chrom" --from-bp "$pos" --to-bp "$end"
gsutil cp "$out_prefix".p{gen,var,sam} ${WORKSPACE_BUCKET}/aryarm/pgens/ALL_SAMPLES/

for EUR_WHITE AFR_BLACK; do
    gsutil cp ${WORKSPACE_BUCKET}/samples/"$pop".csv .
    plink2 --keep <(cut -f1 -d, "$pop".csv | tail -n+2) --out "$out_prefix"."$pop" --pfile "$out_prefix"
    gsutil cp "$out_prefix"."$pop".p{gen,var,sam} ${WORKSPACE_BUCKET}/aryarm/pgens/"$pop"/
done
