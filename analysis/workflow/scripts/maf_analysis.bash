#!/usr/bin/env bash

# arg1: the ID of the causal SNP (ex: 5:88884379)
# arg1: region to extract (ex: 5_87367336-90059999)
# arg2: MAFs (ex: "0.001 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.1")
# arg3: phenotype file (ex: data/aou/phenos/platelet_count.resid.pheno)
# arg4: out path (ex: out)
# ex: workflow/scripts/maf_analysis.bash 5:88884379 5_87367336-90059999 "0.05 0.055 0.1" data/aou/phenos/platelet_count.resid.pheno out
# Note: You should first investigate the MAF of the causal SNP to determine the best MAFs to use:
# plink2 --pfile out/"$region"/genotypes/"$pheno_name"/"$geno_name" --out out/"$region"/genotypes/"$pheno_name"/"$geno_name".maf --freq

set -euo pipefail

best_variant="$1"
region="$2"
maf_str="$3"
read -r -a mafs <<< "$maf_str"
pheno="$4"
out="$5"

output_dir="$out/$region"/mafs
min_maf=0.0001
pheno_name=platelet_count
geno_name=snps.qc.EUR_WHITE
geno_file="$out/$region"/genotypes/"$pheno_name"/"$geno_name"

mkdir -p "$output_dir"

# NOTE: THRESHOLD IS 18 NOT 20
for maf in "${mafs[@]}"; do
    [ ! -f "$output_dir/$maf".hap ] && \
    echo "Running happler for MAF "$maf 1>&2 && \
    happler run \
    -o "$output_dir/$maf".hap \
    --verbosity DEBUG \
    --maf "$maf" \
    --hap-maf "$min_maf" \
    --max-signals 3 \
    --max-iterations 3 \
    --discard-multiallelic \
    --remove-SNPs \
    --indep-thresh 15 \
    -t 18 \
    --chunk-size 500 \
    --out-thresh 5e-08 \
    "$geno_file".pgen \
    "$pheno"
done

echo "Collecting haplotypes with more than one allele" 1>&2
all_mafs=( "${mafs[@]}" )
mafs=$(for maf in "${mafs[@]}"; do echo -e "$(wc -l "$output_dir/$maf".hap)\t$maf"; done | grep -v '^0' | cut -f2)
mafs=( $mafs )

echo "Computing variance explained for each haplotype" 1>&2
workflow/scripts/variance_explained_plot.py \
--verbosity WARNING \
-s <(for maf in "${mafs[@]}"; do echo "$output_dir/$maf".hap; done) \
-o "$output_dir"/variance_explained.png \
"$geno_file".pgen \
"$pheno" \
"$output_dir"/'{maf}'.hap

mkdir -p "$output_dir"/$best_variant
echo "Creating PGEN for best variant" 1>&2
plink2 --snp "$best_variant" --pfile "$geno_file" --make-pgen --freq --out "$output_dir/$best_variant"/best_variant
echo "Transforming all haplotypes into one merged PGEN" 1>&2
# compute LD for each hap at each MAF by merging all of the hap files for each MAF value and transforming them all
haptools transform -o "$output_dir/$best_variant"/haps.pgen "$geno_file".pgen <(
    grep -E '^#' "$output_dir/${mafs[0]}".hap
    for maf in "${mafs[@]}"; do
        sed 's/\tH0\t/\tH0:'"$maf"'\t/;s/\tH1\t/\tH1:'"$maf"'\t/' "$output_dir/$maf".hap | grep -Ev '^#'
    done
)
echo "Computing LD between each haplotype and the causal variant" 1>&2
workflow/scripts/compute_pgen_ld.py --r2 --no-estimate -o "$output_dir/$best_variant"/haps.ld "$output_dir/$best_variant"/haps.pgen "$output_dir/$best_variant"/best_variant.pgen

echo "Computing LD between all SNPs and the causal variant at each threshold" 1>&2
plink2 --r2-unphased 'inter-chr' 'cols=id,freq' --ld-snp "$best_variant" --ld-window-r2 0 --nonfounders --pfile "$geno_file" --out "$output_dir/$best_variant"/snps
echo "Getting the best SNP at each MAF threshold" 1>&2
tail -n+2 "$output_dir/$best_variant"/snps.vcor | sort -k5,5gr | cut -f3-5 > "$output_dir/$best_variant"/snps.sort.vcor
set +o pipefail
# create an r^2 report for the SNPs
{
    echo -e "maf_thresh\tsnp\tmaf\tr2"
    for maf in "${mafs[@]}"; do
        echo -ne "$maf\t"
        awk -F $'\t' '$2 > '"$maf" "$output_dir/$best_variant"/snps.sort.vcor | head -1
    done
} > "$output_dir/$best_variant"/snps.ld
set -o pipefail

echo "Merging the SNP and hap r2 reports together" 1>&2
{
    echo -e "maf_thresh\thap_id\thap_r2\tsnp_r2"
    join -t $'\t' -12 -21 <(
        tail -n+2 "$output_dir/$best_variant"/haps.ld | cut -f3,4 | sed 's/:/\t/' | sort -k2,2
    ) <(
        tail -n+2 "$output_dir/$best_variant"/snps.ld | cut -f 1,4 | sort -k1,1
    )
} > "$output_dir/$best_variant"/maf_hap_snp_r2.tsv

# -------------

cd "$output_dir"

echo "Plotting variance explained" 1>&2
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

cd "$best_variant"

echo "Plotting LD with causal variant" 1>&2
(
  echo "a=["$(tail -n+2  maf_hap_snp_r2.tsv | sed 's/\tH0\t/\t0\t/;s/\tH1\t/\t1\t/' | sort -t$'\t' -k1,1g | tr $'\t' , | sed 's/^/(/;s/$/)/' | paste -s -d,)"]"
  cat <<'EOF'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = np.array(a)
hap_ids = np.unique(data[:,1])
for hp_id in hap_ids:
  dat = data[hp_id == data[:,1]]
  plt.plot(-np.log10(dat[:,0]), dat[:,2], 'o', label=f"Haplotype {int(hp_id)}")
plt.plot(-np.log10(data[:,0]), data[:,3], 'o', label="Best SNP")
plt.xlabel("-log10(MAF)")
plt.ylabel("LD (R^2) with best SNP at --mac 20")
plt.legend()
plt.savefig("ld_maf.png")

EOF
) | python
