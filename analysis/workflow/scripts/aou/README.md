# scripts
This directory contains various scripts used to prepare data in AoU.

Any python scripts implement the --help argument. For bash, awk, and R scripts, you can run `head <script>` to read about their usage.

## [get_pgen.bash](get_pgen.bash)
This script generates a PGEN file for a specific locus by subsetting and merging a series of phased VCFs in our AoU v7 bucket. It is used within the aou_v7 rule within the Snakemake workflow.

## [get_phenos.bash](get_phenos.bash)
TODO

## [merge_batched_vcfs.bash](merge_batched_vcfs.bash)
TODO

## [plink2_qc_EUR_AFR.bash](plink2_qc_EUR_AFR.bash)
TODO

## [vcf2pgen.bash](vcf2pgen.bash)
TODO
