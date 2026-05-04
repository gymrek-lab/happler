#!/usr/bin/env bash

# arg1: the ID of the causal SNP (ex: 5:88884379)
# arg1: region to extract (ex: 5_87367336-90059999)
# arg2: MAFs (ex: "0.001 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.065 0.085 0.1")
# arg3: phenotype file (ex: data/aou/phenos/platelet_count.resid.pheno)
# arg4: out path (ex: out)
# ex: workflow/scripts/maf_analysis.bash 5:88884379 5_87367336-90059999 "0.05 0.055 0.1" data/aou/phenos/platelet_count.resid.pheno out
# Note: You should first investigate the MAF of the causal SNP to determine the best MAFs to use:
# plink2 --pfile out/"$region"/genotypes/"$pheno_name"/"$geno_name" --out out/"$region"/genotypes/"$pheno_name"/"$geno_name".maf --freq

eval "$(conda shell.bash hook)"
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
    "$pheno" &> "$output_dir/$maf".log
done

echo "Collecting haplotypes with more than one allele" 1>&2
all_mafs=( "${mafs[@]}" )
mafs=$(for maf in "${mafs[@]}"; do echo -e "$(wc -l "$output_dir/$maf".hap)\t$maf"; done | grep -v '^0' | cut -f2)
mafs=( $mafs )

echo "Transforming all haplotypes into one merged PGEN" 1>&2
# merge all of the hap files for each MAF value and transform them all
haptools transform -o "$output_dir"/haps.pgen "$geno_file".pgen <(
    grep -E '^#' "$output_dir/${mafs[0]}".hap
    for maf in "${mafs[@]}"; do
        sed 's/\tH0\t/\tH0:'"$maf"'\t/;s/\tH1\t/\tH1:'"$maf"'\t/' "$output_dir/$maf".hap | grep -Ev '^#'
    done
) &> "$output_dir"/haps.log

for maf in "${all_mafs[@]}"; do
    mkdir -p "$output_dir/susie_$maf"
    if printf '%s\0' "${mafs[@]}" | grep -qzxF "$maf"; then
        echo "Merging SNPs and haps for MAF $maf" 1>&2
        workflow/scripts/merge_plink.py \
        --chunk-size 1000 \
        --maf "$maf" \
        --maf-file 2 \
        --extract <(grep -Ev '^#' "$output_dir"/haps.pvar | cut -f3 | grep ':'"$maf") \
        --verbosity DEBUG \
        "$output_dir"/haps.pgen \
        "$geno_file".pgen \
        "$output_dir/susie_$maf"/merge.pgen &> "$output_dir/susie_$maf"/merge.log
    else
        echo "Filtering SNPs for MAF $maf" 1>&2
        plink2 --maf "$maf" --pfile "$geno_file" --make-pgen --out "$output_dir/susie_$maf"/merge &>/dev/null
    fi
    conda activate .snakemake/conda/885b27680699bdbf0ec4008de1a842a3_
    [ ! -f "$output_dir/susie_$maf"/susie.rds ] && \
    echo "Running SuSiE for MAF $maf" 1>&2 && \
    workflow/scripts/run_SuSiE.R "$output_dir/susie_$maf"/merge.pgen "$pheno" "$output_dir/susie_$maf" NULL "$(echo "$region" | sed 's/_/:/')" 10 &> "$output_dir/susie_$maf"/susie.log && \
    workflow/scripts/extract_pips.R "$output_dir/susie_$maf"/susie.rds "$output_dir/susie_$maf"/pips.tsv &>"$output_dir/susie_$maf"/pips.log
    conda deactivate
done

echo "Computing variance explained for each haplotype" 1>&2
workflow/scripts/variance_explained_plot.py \
--verbosity WARNING \
-s <(for maf in "${mafs[@]}"; do echo "$output_dir/$maf".hap; done) \
-o "$output_dir"/variance_explained.png \
"$geno_file".pgen \
"$pheno" \
"$output_dir"/'{maf}'.hap &> "$output_dir"/variance_explained.log

mkdir -p "$output_dir"/$best_variant
echo "Creating PGEN for best variant" 1>&2
plink2 --snp "$best_variant" --pfile "$geno_file" --make-pgen --freq --out "$output_dir/$best_variant"/best_variant &>/dev/null
best_variant_maf="$(grep -P "\t$best_variant\t" "$output_dir/$best_variant"/best_variant.afreq | cut -f5 | awk '{min = ($1 < 1-$1 ? $1 : 1-$1); print min;}')"
echo "Causal variant MAF: $best_variant_maf" 1>&2
echo "Computing LD between each haplotype and the causal variant" 1>&2
workflow/scripts/compute_pgen_ld.py --r2 --no-estimate -o "$output_dir/$best_variant"/haps.ld "$output_dir"/haps.pgen "$output_dir/$best_variant"/best_variant.pgen &> "$output_dir/$best_variant"/haps.ld.log

echo "Computing LD between all SNPs and the causal variant at each threshold" 1>&2
plink2 --r2-unphased 'inter-chr' 'cols=id,freq' --ld-snp "$best_variant" --ld-window-r2 0 --nonfounders --pfile "$geno_file" --out "$output_dir/$best_variant"/snps
echo "Getting the best SNP at each MAF threshold" 1>&2
tail -n+2 "$output_dir/$best_variant"/snps.vcor | sort -k5,5gr | cut -f3-5 > "$output_dir/$best_variant"/snps.sort.vcor

set +o pipefail
# create an r^2 report for the SNPs
{
    echo -e "maf_thresh\tsnp\tmaf\tr2"
    for maf in "${all_mafs[@]}"; do
        echo -ne "$maf\t"
        if awk -v n1="$maf" -v n2="$best_variant_maf" 'BEGIN { exit (n1 > n2 ? 0 : 1) }'; then
            awk -F $'\t' '$2 > '"$maf" "$output_dir/$best_variant"/snps.sort.vcor | head -1
        else
            echo -e "$best_variant\t$best_variant_maf\t1"
        fi
    done
} > "$output_dir/$best_variant"/snps.ld
set -o pipefail

echo "Merging the SNP and hap r2 reports together" 1>&2
{
    echo -e "maf_thresh\thap_id\thap_r2\tsnp_r2\tsnp_maf"
    join -t $'\t' -a1 -11 -22 <(
        tail -n+2 "$output_dir/$best_variant"/snps.ld | cut -f 1,3,4 | sort -k1,1
    ) <(
        tail -n+2 "$output_dir/$best_variant"/haps.ld | cut -f3,4 | sed 's/:/\t/' | sort -k2,2
    ) | \
    awk -F $'\t' -v 'OFS=\t' 'NF == 3 { $4="H0";$5=$3 } {print $1, $4, $5, $3, $2;}'
} > "$output_dir/$best_variant"/maf_hap_snp_r2.tsv

# -------------

echo "Plotting variance explained" 1>&2
( cd "$output_dir" && (
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
) | python; )

echo "Plotting LD with causal variant" 1>&2
( cd "$output_dir/$best_variant" && (
  echo "causal_variant=\"$best_variant\"; a=["$(tail -n+2 maf_hap_snp_r2.tsv | sed 's/\tH0\t/\t0\t/;s/\tH1\t/\t1\t/' | sort -t$'\t' -k1,1g | tr $'\t' , | sed 's/^/(/;s/$/)/' | paste -s -d,)"]"
  cat <<'EOF'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = np.array(a)
both = data[data[:,2] == data[:,3]]
hap_ids = np.unique(data[:,1])
lowest_maf = np.unique(data[:,4]).min()
# Pick a colormap and sample as many distinct colors as needed
cmap = plt.get_cmap("tab20")  # good for up to ~20 distinct colors
colors = cmap(np.linspace(0, 1, len(hap_ids)))
color_by_hap = {hp: colors[i] for i, hp in enumerate(hap_ids)}
for hp_id in hap_ids:
    dat = data[hp_id == data[:,1]]
    plt.plot(dat[:,0], dat[:,2], 'o-', color=color_by_hap[hp_id], label=f"Haplotype {int(hp_id)}")
plt.axvline(x=lowest_maf, color='red', linestyle='--')
plt.plot(data[:,0], data[:,3], 'o-', color='black', label="Best SNP")
plt.plot(both[:,0], both[:,3], 'o', color='grey', label="Both")
plt.xlabel("MAF")
plt.ylabel("LD (R^2) with causal SNP")
plt.ylim(0, 1.02)
plt.title(causal_variant)
plt.legend()
plt.savefig("ld_maf.png")

EOF
) | python; )
