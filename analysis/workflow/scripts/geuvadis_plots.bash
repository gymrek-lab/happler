#!/usr/bin/env bash
# Recreates a few useful plots for interpreting the output of the Geuvadis analysis
# arg1: The path to the "out" folder

out="$1"





############################################## MAIN PROGRAM ########################################

# first, create the multiline.txt file, which lists all .hap files with substantial haplotypes
while read hap; do ls "$out/$(echo "$hap" | cut -f1)/happler/run/$(echo "$hap" | cut -f2)/happler.hap"; done < $out/multiline.tsv > $out/multiline.txt

# let's report a few statistics
num_tot_regions="$(ls -d $out/*_*-* | wc -l)"
num_regions="$(cat $out/multiline.txt | wc -l)"
echo "Out of $num_tot_regions regions, $num_regions has at least one haplotype with more than one variant."
echo "Of those $num_regions, here is a breakdown of the number of haplotypes each region had:"
while read hap; do grep '^H' $hap | wc -l; done < $out/multiline.txt | sort | uniq -c
avg_num_alleles="$(for i in $(cat $out/multiline.txt); do grep '^V' $i | cut -f2 | sort | uniq -c | sed 's/^ *//' | cut -f1 -d' '; done | awk '{ total += $1 } END { print total/NR }')"
echo "Of those $num_regions, the average number of alleles in each haplotype is $avg_num_alleles."

# now, let's make the variance_explained.png plot
workflow/scripts/variance_explained_plot.py --verbosity WARNING -o "$out/variance_explained.png" -s "$out"/multiline.txt "$out"/{locus}/happler/run/{gene}/include/merged.pgen data/geuvadis/phenos/{gene}.pheno "$out"/{locus}/happler/run/{gene}/happler.hap
echo "Created $out/variance_explained.png"

# now, let's make the pips.tsv file
for i in "$out"/*/happler/run/*/include/susie_pips.tsv; do grep -P '^H[0-9]\t' $i | sed 's+^+'"$(echo $i | sed 's\/include/susie_pips.tsv$\\;s\^out/\\;s+/happler/run/+:+')"':+'; done > "$out"/pips.tsv
echo "Created $out/pips.tsv"
echo "$(awk -F $'\t' '$2 > 0.9' "$out"/pips.tsv | wc -l) of the $(cat "$out"/pips.tsv | wc -l) haplotypes have PIPs above 0.9"

cd "$out"

# now, let's make the hap_pips.png file
(
  echo "a=["$(cut -f2 pips.tsv | paste -s -d,)"]"
  cat <<'EOF'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = np.array(a)
binwidth = 0.1
plt.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth))
plt.title("Haplotype PIPs")
plt.tight_layout()
plt.savefig("hap_pips.png", bbox_inches="tight")
EOF
) | python
echo "Created $out/hap_pips.png"

# now, let's check the HWE of the haplotypes
for i in $(cat multiline.txt | sed 's/.hap$/.pvar/;s+out/++'); do plink2 --pfile ${i%.*} --hardy --out ${i%.*}-hwe &>/dev/null; done
(
  echo "a=["$(for i in */happler/run/*/happler-hwe.hardy; do cut -f10 $i | tail -n+2; done | paste -s -d,)"]"
  cat <<'EOF'
import numpy as np
import matplotlib.pyplot as plt
data = -10*np.log10(np.array(a))
binwidth = 2
plt.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth))
plt.title("Haplotype HWE -10log10 P-Values")
plt.savefig("hwe.png")
EOF
) | python
echo "Created $out/hwe.png"

# now, let's threshold by MAC and create a histogram
for i in $(cat multiline.txt | sed 's/.hap$/.pvar/;s+out/++'); do plink2 --pfile ${i%.*} --mac 70 --make-pgen --out ${i%.*}-maf --freq &>/dev/null; done
echo "$(grep 'Error: No variants remaining' */happler/run/*/happler-maf.log | cut -d '/' -f2,5 | wc -l) haplotypes had an MAC below 70."

# create SV LD plot
# first, copy all of the results over
mkdir -p sv_ld/H0
for i in */happler/run/*/happler_svs.ld; do region="$(echo "$i" | sed 's\/happler_svs.ld$\\;s+/happler/run/+\t+')"; cp "$i" sv_ld/H0/$(echo "$region" | cut -f1):$(echo "$region" | cut -f2).ld; done
# now, collate the results
{ echo -e 'file\tpip\tpos\tid\tld'; sort -gr -k2,2 pips.tsv | { while read -r line; do file="sv_ld/$(echo "$line" | cut -f1 | cut -f3 -d:)/$(echo "$line" | cut -f1 | cut -f-2 -d:).ld"; echo -en "$file"$'\t'; echo -en "$(echo "$line" | cut -f2)"$'\t'; awk -F $'\t' -v 'OFS=\t' '{print $2, $3, sqrt($4*$4);}' "$file" | sort -gr -k3,3 | head -n1; done } | sed 's/^.*sv_ld\///'; } > pips_sv_ld.tsv
echo "Created $out/pips_sv_ld.tsv"
# now, visualize all of the results
(
  echo "a=["$(tail -n+2 pips_sv_ld.tsv | cut -f2,5 | tr $'\t' , | sed 's/^/(/;s/$/)/' | paste -s -d,)"]"
  cat <<'EOF'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = np.array(a)
plt.scatter(data[:,0], data[:,1])
plt.xlim(-0.02, 1.02)
plt.ylim(-0.02, 1.02)
plt.axline([0, 0], [1, 1])
plt.gca().xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
plt.gca().yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
plt.grid(True, which='both')
plt.xlabel("Haplotype PIP")
plt.ylabel("LD of Best SV")
plt.savefig("pips_best_sv_ld.png")

EOF
) | python
echo "Created $out/pips_best_sv_ld.png"

cd -
# create STR LD plot
# first, let's get a list of the STRs in the regions with haplotypes
(
  echo -ne "region\tgene\thap.id\tpip\t"
  head -n1 data/geuvadis/mlamkin/Geuvadis_varlevel_corrected_significant_variants.with-end.tsv && \
  ~/miniconda3/envs/htslib/bin/bedtools intersect -a <(
    echo -e "chrom\tstart\tend\tgene\tpip" && cat "$out"/pips.tsv | sed 's+_+\t+;s+-+\t+;s+:+\t+g' | sort -k1,1V -k2,2n
  ) -b data/geuvadis/mlamkin/Geuvadis_varlevel_corrected_significant_variants.with-end.tsv -wa -wb -loj | \
  sed 's+\t+_+;s+\t+-+'
) | awk -F'\t' '$2 == $8' | cut -f8 --complement > "$out"/STR_assocations.tsv
echo "Created $out/STR_assocations.tsv"
# now, let's compute LD for each region
echo -e "hap\tpip\tpos\tid\tld\talleles" > "$out"/pips_str_ld.tsv
while IFS= read -r line; do
  str_id="$(echo "$line" | cut -f5,6 --output-delimiter ':')"
  echo -ne "$(echo "$line" | cut -f1-3 --output-delimiter ':')\t$(echo "$line" | cut -f4)\t$(echo "$line" | cut -f6)\t$str_id\t"
  echo -ne "$(workflow/scripts/compute_pgen_ld.py --verbosity WARNING --target-is-repeat --hap-id "$str_id" -o /dev/stdout "$out/$(echo "$line" | cut -f1)"/happler/run/"$(echo "$line" | cut -f2)"/happler.pgen data/geuvadis/mlamkin/all_Geuvadis_STRs.pgen | tail -n+2 | cut -f4)"
  echo -e "\t$(grep -P '\t'"$str_id"'\t' data/geuvadis/mlamkin/all_Geuvadis_STRs.pvar | cut -f 4,5 --output-delimiter ,)"
done < <(tail -n+2 "$out"/STR_assocations.tsv) >> "$out"/pips_str_ld.tsv
echo "Created $out/pips_str_ld.tsv"
(
  head -n1 "$out/pips_str_ld.tsv"
  tail -n+2 "$out/pips_str_ld.tsv" \
    | awk -F'\t' -v OFS='\t' '{$5 = ($5 < 0) ? -$5 : $5; print}' \
    | sort -t$'\t' -k1,1 -k5,5nr \
    | awk -F'\t' -v OFS='\t' '!seen[$1]++ { print }'
) > "$out/pips_best_str_ld.tsv"
echo "Created $out/pips_best_str_ld.tsv"
cd "$out"
(
  echo "a=["$(tail -n+2 pips_best_str_ld.tsv | cut -f2,5 | tr $'\t' , | sed 's/^/(/;s/$/)/' | paste -s -d,)"]"
  cat <<'EOF'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = np.array(a)
plt.scatter(data[:,0], data[:,1])
plt.xlim(-0.02, 1.02)
plt.ylim(-0.02, 1.02)
plt.axline([0, 0], [1, 1])
plt.gca().xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
plt.gca().yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
plt.grid(True, which='both')
plt.xlabel("Haplotype PIP")
plt.ylabel("LD of Best STR")
plt.savefig("pips_best_str_ld.png")

EOF
) | python
echo "Created $out/pips_best_str_ld.png"
# now, compare STR vs SV LD
(
  echo "a=["$(join -t $'\t' -j1 --header <(head -n1 pips_sv_ld.tsv; tail -n+2 pips_sv_ld.tsv | sed 's/.ld\t/:H0\t/;s+^H0/++' | sort -k1,1) <(head -n1 pips_best_str_ld.tsv; tail -n+2 pips_best_str_ld.tsv | sort -k1,1) | cut -f5,9 | tail -n+2 | tr $'\t' , | sed 's/^/(/;s/$/)/' | paste -s -d,)"]"
  cat <<'EOF'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = np.array(a)
# plt.scatter(x, y)
plt.scatter(data[:,0], data[:,1])
plt.xlim(-0.02, 1.02)
plt.ylim(-0.02, 1.02)
plt.axline([0, 0], [1, 1])
plt.gca().xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
plt.gca().yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
plt.grid(True, which='both')
plt.xlabel("LD of best SV")
plt.ylabel("LD of Best STR")
plt.savefig("ld_str_vs_sv.png")

EOF
) | python
echo "Created $out/ld_str_vs_sv.png"
