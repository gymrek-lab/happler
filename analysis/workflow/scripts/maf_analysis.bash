#!/usr/bin/env bash

# arg1: the ID of the causal SNP (ex: 5:88884379)
# arg1: region to extract (ex: 5_87367336-90059999)
# arg2: MAFs (ex: (0.001 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.1))
# arg3: phenotype file (ex: data/aou/phenos/platelet_count.resid.pheno)
# arg4: out path (ex: out)
# ex: workflow/scripts/maf_analysis.bash 5:88884379 5_87367336-90059999 (0.05 0.055 0.1) data/aou/phenos/platelet_count.resid.pheno out

# Note: You should first investigate the MAF of the causal SNP to determine the best MAFs to use:
# plink2 --pfile out/"$region"/genotypes/"$pheno_name"/"$geno_name" --out out/"$region"/genotypes/"$pheno_name"/"$geno_name".maf --freq

best_variant="$1"
region="$2"
mafs=("${3[@]}")
pheno="$4"
out="$5"

output_dir="$out/$region"/mafs
min_maf=0.0001
pheno_name=platelet_count
geno_name=snps.qc.EUR_WHITE

mkdir -p "$output_dir"

for maf in "${mafs[@]}"; do
    happler run \
    -o "$output_dir"/"$maf".hap \
    --verbosity DEBUG \
    --maf "$maf" \
    --hap-maf "$min_maf" \
    --max-signals 3 \
    --max-iterations 3 \
    --discard-multiallelic \
    --remove-SNPs \
    --indep-thresh 15 \
    -t 20 \
    --chunk-size 500 \
    --out-thresh 5e-08 \
    "$out/$region"/genotypes/"$pheno_name"/"$geno_name".pgen \
    "$pheno"
done

workflow/scripts/variance_explained_plot.py \
--verbosity WARNING \
-s <(wc -l "$output_dir"/*.hap | grep -v total | sed 's/^ *//' | grep -v '^0' | sed 's/^.* //') \
-o "$output_dir"/variance_explained.png \
"$out/$region"/genotypes/"$pheno_name"/"$geno_name".pgen \
"$pheno" \
"$output_dir"/{maf}.hap

# compute LD for each hap at each MAF by merging all of the hap files for each MAF value and transforming them all
haptools transform -o "$output_dir"/$best_variant/haps.pgen "$out/$region"/genotypes/"$pheno_name"/"$geno_name".pgen <(grep -E '^#' "$output_dir"/"${mafs[0]}".hap; for maf in "${mafs[@]}"; do sed 's/\tH0\t/\tH0:'"$maf"'\t/;s/\tH1\t/\tH1:'"$maf"'\t/' "$output_dir"/$maf.hap | grep -Ev '^#'; done)
workflow/scripts/compute_pgen_ld.py --r2 --no-estimate -o "$output_dir"/$best_variant/haps.ld "$output_dir"/$best_variant/haps.pgen "$output_dir"/$best_variant/best_variant.pgen

cd "$output_dir"

# plot variance_explained
(
  echo "a=["$(cut -f1,6 variance_explained.tsv | tail -n+2 | tr ':' $'\t' | grep 'H0' | cut -f 1,3 | sort -g | tr $'\t' , | sed 's/^/(/;s/$/)/' | paste -s -d,)"]"
  cat <<'EOF'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = np.array(a)
plt.plot(-np.log10(data[:,0]), data[:,1], '-o')
plt.xlabel("-log10(MAF)")
plt.ylabel("Haplotype / Haplotype's SNPs")
plt.title("Variance Explained (R^2)")
plt.savefig("varexp_maf.png")

EOF
) | python
