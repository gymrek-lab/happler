#!/usr/bin/env bash

# Convert the phased AoU v8 VCFs to PGENs
# Allocate 4 CPUs, 26 GB memory, and a 400 GB SSD (total: $10.19 for 24 hrs)
# ETA: 7-10 days

V8_BUCKET="$1" # ex: gs://fc-secure-MYBUCKET
V8_CDR_DIR="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/phasing"

mkdir -p ~/phased_pgens
cd ~/phased_pgens

for i in {22..1}; do
    gsutil -u $GOOGLE_PROJECT -m cp -r "$V8_CDR_DIR"/chr${i}_*.vcf.gz chr${i}.vcf.gz && \
    plink2 --memory 24000 --vcf chr${i}.vcf.gz --out chr${i} --set-all-var-ids '@:#' --make-just-pvar --keep-autoconv && \
    rm chr${i}.vcf.gz;
    (
        gsutil cp chr${i}.log chr${i}.p* "$V8_BUCKET"/phased_pgens/ && \
        rm chr${i}.log chr${i}.p* && \
        echo chr${i};
    ) &
done
wait
echo finally
