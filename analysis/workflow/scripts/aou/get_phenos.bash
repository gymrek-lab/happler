#!/usr/bin/env bash

export GCS_REQUESTER_PAYS_PROJECT="${GOOGLE_PROJECT}"
export GCS_OAUTH_TOKEN="$(gcloud auth application-default print-access-token)"
CDR_DIR="gs://fc-aou-datasets-controlled/v8"

mkdir -p phenos
cd phenos

gsutil -u "$GOOGLE_PROJECT" cat "$CDR_DIR/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv" | \
cut -f1,4 | \
sed 's/\[/\t/;s/\]/\t/;s/, /\t/g' | \
tail -n+2 | {
    echo -ne "#IID\t\t"
    read -r head
    seq 1 "$(echo "$head" | cut -f 3- | awk -F $'\t' '{print NF; exit}')" | \
    sed 's/^/PC/' | \
    paste -s -d$'\t'; echo "$head"; cut --complement -f17
} | \
cut --complement -f2 > AOU_PCS.covar

for pheno in platelet_count_new_phenocovar ldl_cholesterol_phenocovar; do
    gsutil cp ${WORKSPACE_BUCKET}/phenotypes/$pheno.csv $pheno.csv
    cut -f-2 -d, --output-delimiter $'\t' $pheno.csv | { echo -e "#IID\tpheno"; tail -n+2; } > $pheno.og.pheno
    cut -f 1,3- -d, --output-delimiter $'\t' $pheno.csv | { read -r head; echo "$head" | sed 's/person_id/#IID/'; cat; } > $pheno.og.covar
    # we use the first 10 PCs and the first column is #IID, so we do -f-11
    ../happler/analysis/workflow/scripts/residuals.py -o $pheno.residuals.pheno -e <(cut -f-11 AOU_PCS.covar) $pheno.og.pheno $pheno.og.covar
    gsutil cp $pheno.residuals.pheno ${WORKSPACE_BUCKET}/aryarm/data/aou/phenos/
done
