#!/usr/bin/env bash

region="$1"
pop="${2:-EUR_WHITE}"

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
    bcftools view -O z -o "$batch".bcf -r "$region" "$(echo "$VCF_DIR" | sed 's/*/'"$batch"'/')"
#    plink2 --out "$batch" --nonfounders --bcf "$batch".bcf --geno 0 --make-pgen --allow-extra-chr --max-alleles 2 --chr "$chrom" --from-bp "$pos" --to-bp "$end"
done

cd ..   
bcftools merge --no-index -O z -o "$out_prefix".bcf -l <(ls "$out_prefix"/*.bcf)
# note that we skip --maf bc the input is already filtered
plink2 --out "$out_prefix" --nonfounders --bcf "$out_prefix".bcf --geno 0 --make-pgen --allow-extra-chr --max-alleles 2 --chr "$chrom" --from-bp "$pos" --to-bp "$end"
gsutil cp "$out_prefix".p{gen,var,sam} ${WORKSPACE_BUCKET}/aryarm/pgens/ALL_SAMPLES/

gsutil cp ${WORKSPACE_BUCKET}/samples/"$pop".csv .
plink2 --keep <(cut -f1 -d, "$pop".csv | tail -n+2) --out "$out_prefix"."$pop" --pfile "$out_prefix"
gsutil cp "$out_prefix"."$pop".p{gen,var,sam} ${WORKSPACE_BUCKET}/aryarm/pgens/"$pop"/

exit
# to process phenotypes as well
mkdir -p phenos
gsutil -u "$GOOGLE_PROJECT" cat "$CDR_DIR/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv" | cut -f1,4 |  sed 's/\[/\t/;s/\]/\t/;s/, /\t/g' | tail -n+2 | { echo -ne "#IID\t\t"; read -r head; seq 1 "$(echo "$head" | cut -f 3- | awk -F $'\t' '{print NF; exit}')" | sed 's/^/PC/' | paste -s -d$'\t'; echo "$head"; cut --complement -f17; } | cut --complement -f2 > phenos/AOU_PCS.covar
for pheno in platelet_count_new_phenocovar ldl_cholesterol_phenocovar; do
    gsutil cp ${WORKSPACE_BUCKET}/phenotypes/$pheno.csv phenos/$pheno.csv
    cut -f-2 -d, --output-delimiter $'\t' phenos/$pheno.csv | { echo -e "#IID\tpheno"; tail -n+2; } > phenos/$pheno.og.pheno
    cut -f 1,3- -d, --output-delimiter $'\t' phenos/$pheno.csv | { read -r head; echo "$head" | sed 's/person_id/#IID/'; cat; } > phenos/$pheno.og.covar
    happler/analysis/workflow/scripts/residuals.py -o phenos/$pheno.residuals.pheno -e <(cut -f-11 phenos/AOU_PCS.covar) phenos/$pheno.og.pheno phenos/$pheno.og.covar
done
