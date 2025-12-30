#!/usr/bin/env bash

# A Bash script using plink2 to perform genotype QC with EUR, AFR, and phenotype-defined cohorts,
# while ensuring multiallelic variant filtering is applied at the final filtering step.

set -euo pipefail

# USAGE: hail_qc_EUR_AFR.sh <genotypes.pgen> <phenotypes.pheno> [-o output.pgen]
# Optional arguments:
#   --samples-file-dir       Path to EUR_WHITE.csv and AFR_BLACK.csv files (default: current directory)
#   --sample-call-rate       Minimum sample call rate (default: 0.9)
#   --variant-call-rate      Minimum variant call rate (default: 0.9)
#   --MAF                    Minor allele frequency QC threshold (default: 0.01)
#   --HWE                    HWE p-value cutoff QC threshold (default: 1e-100)

# Default values
SAMPLES_DIR="."
SAMPLE_CALL_RATE=0.9
VARIANT_CALL_RATE=0.9
MAF=0.01
HWE=1e-100

# Parse arguments
GENOTYPES="$1"
PHENOTYPES="$2"
shift 2

OUTPUT=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        --samples-file-dir)
            SAMPLES_DIR="$2"
            shift 2
            ;;
        --sample-call-rate)
            SAMPLE_CALL_RATE="$2"
            shift 2
            ;;
        --variant-call-rate)
            VARIANT_CALL_RATE="$2"
            shift 2
            ;;
        --MAF)
            MAF="$2"
            shift 2
            ;;
        --HWE)
            HWE="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [ -z "$OUTPUT" ]; then
    OUTPUT="${GENOTYPES%.pgen}.qc"
fi

# Ensure required inputs exist
test -f "$GENOTYPES" || { echo "Genotypes file $GENOTYPES not found"; exit 1; }
test -f "$PHENOTYPES" || { echo "Phenotype file $PHENOTYPES not found"; exit 1; }

# Ensure the presence of sample cohort files
EUR_FILE="${SAMPLES_DIR}/EUR_WHITE.csv"
AFR_FILE="${SAMPLES_DIR}/AFR_BLACK.csv"
test -f "$EUR_FILE" || { echo "EUR cohort file ${EUR_FILE} not found"; exit 1; }
test -f "$AFR_FILE" || { echo "AFR cohort file ${AFR_FILE} not found"; exit 1; }

# Convert cohort files to space-delimited sample lists
cut -d, -f1 "$EUR_FILE" | tail -n +2 > eur_samples.keep
cut -d, -f1 "$AFR_FILE" | tail -n +2 > afr_samples.keep

# Extract phenotype samples as its own cohort
cut -f1 "$PHENOTYPES" | tail -n +2 > pheno_samples.keep

# Step 1: Compute HWE and AF in the EUR cohort
plink2 --pfile "${GENOTYPES%.pgen}" --keep eur_samples.keep \
       --freq --hardy --out eur_stats

# Extract EUR variants passing the criteria
awk -v maf="$MAF" -v hwe="$HWE" '$5 >= maf && $8 > hwe {print $2}' eur_stats.hardy > eur_passing_variants.txt

# Step 2: Compute HWE and AF in the AFR cohort
plink2 --pfile "${GENOTYPES%.pgen}" --keep afr_samples.keep \
       --freq --hardy --out afr_stats

# Extract AFR variants passing the criteria
awk -v maf="$MAF" -v hwe="$HWE" '$5 >= maf && $8 > hwe {print $2}' afr_stats.hardy > afr_passing_variants.txt

# Step 3: Compute HWE and AF in the phenotype-defined cohort
plink2 --pfile "${GENOTYPES%.pgen}" --keep pheno_samples.keep \
       --freq --hardy --out pheno_stats

# Extract PHENO variants passing the criteria
awk -v maf="$MAF" -v hwe="$HWE" '$5 >= maf && $8 > hwe {print $2}' pheno_stats.hardy > pheno_passing_variants.txt

# Step 4: Combine variant lists (union of EUR, AFR, and PHENO passing variants)
cat eur_passing_variants.txt afr_passing_variants.txt pheno_passing_variants.txt | sort | uniq > passing_variants.txt

# Step 5: Filter the original pgen file using the list of passing variants and remove multiallelics
plink2 --pfile "${GENOTYPES%.pgen}" --extract passing_variants.txt \
       --max-alleles 2 --mind "$SAMPLE_CALL_RATE" --geno "$VARIANT_CALL_RATE" \
       --make-pgen --out "$OUTPUT"

# Final output
echo "Filtered genotype data saved to ${OUTPUT}.pgen"
