#!/usr/bin/env bash

# Regress PCs and other covariates out of AoU phenotypes
# Allocate 4 CPUs, 3.6 GB memory, and a 120 GB disk (the cheapest possible configuration)
# Execute this script from within the home directory

CDR_DIR="gs://fc-aou-datasets-controlled/v8"

mkdir -p phenos
cd phenos

gcloud storage --billing-project "$GOOGLE_PROJECT" cat "$CDR_DIR/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv" | \
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
    gcloud storage cp ${WORKSPACE_BUCKET}/phenotypes/$pheno.csv $pheno.csv
    cut -f-2 -d, --output-delimiter $'\t' $pheno.csv | { echo -e "#IID\tpheno"; tail -n+2; } > $pheno.og.pheno
    cut -f 1,3- -d, --output-delimiter $'\t' $pheno.csv | { read -r head; echo "$head" | sed 's/person_id/#IID/'; cat; } > $pheno.og.covar
    # we use the first 10 PCs and the first column is #IID, so we do -f-11
    ../happler/analysis/workflow/scripts/residuals.py -o $pheno.residuals.pheno -e <(cut -f-11 AOU_PCS.covar) $pheno.og.pheno $pheno.og.covar
    # Note that we passed both .covar files. If we wanted to merge the .covar files instead, we could do it like this:
    # join -j 1 -t $'\t' --header <(cat $pheno.og.covar | (sed-u 1q; sort -k1,1)) <(cat AOU_PCS.covar | cut -f-11 | (sed -u 1q; sort -k1,1)) > $pheno.merged.og.covar
    gcloud storage cp $pheno.residuals.pheno ${WORKSPACE_BUCKET}/aryarm/data/aou/phenos/$(echo $pheno | sed 's/_new_/_/').resid.pheno
done
