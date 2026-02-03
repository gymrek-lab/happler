# scripts
This directory contains various scripts used to prepare data in AoU.

Any python scripts implement the --help argument. For bash, awk, and R scripts, you can run `head <script>` to read about their usage.

## [get_pgen.bash](get_pgen.bash)
Generate a PGEN file for a specific locus by subsetting and merging a series of phased VCFs in our AoU v7 bucket. This script is used within the aou_v7 rule within the Snakemake workflow.

## [get_phenos.bash](get_phenos.bash)
Regress PCs and other covariates out of AoU phenotypes and upload the result to our storage bucket.

## [merge_batched_vcfs.bash](merge_batched_vcfs.bash)
Efficiently merge VCFs which share variants but have non-overlapping sample sets using paste and not bcftools

## [plink2_qc_EUR_AFR.bash](plink2_qc_EUR_AFR.bash)
Perform QC steps on AoU PGENs. This script was used within the Snakemake workflow but was replaced by the genotypes_aou_qc rule and is now no longer used.

## [vcf2pgen.bash](vcf2pgen.bash)
Convert the phased AoU v8 VCFs to PGENs and reupload them to our storage bucket.
