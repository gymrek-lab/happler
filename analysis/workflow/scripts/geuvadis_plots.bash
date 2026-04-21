#!/usr/bin/env bash
# Recreates a few useful plots for interpreting the output of the Geuvadis analysis
# arg1: The path to the "out" folder
# arg2: The operating "mode" (either "geuvadis", "ukb", or "aou")

# Required inputs:
# out/{locus}/happler/run/{gene}/bench/run
# out/{locus}/happler/run/{gene}/happler.hap
# out/{locus}/happler/run/{gene}/happler.p{gen,var,sam}
# out/{locus}/happler/run/{gene}/include/susie_pips.tsv
# out/{locus}/happler/run/{gene}/exclude/susie_pips.tsv
# out/{locus}/happler/run/{gene}/include/merged.pgen (could be switched out for the original snp panels ie snps.pgen)
# out/{locus}/genotypes/{gene}/snps.pgen

out="${1:out}"
mode="${2:aou}"




############################################## MAIN PROGRAM ########################################

# first, create the multiline.txt file, which lists all .hap files with substantial haplotypes
while read hap; do echo "out/$(echo "$hap" | cut -f1)/happler/run/$(echo "$hap" | cut -f2)/happler.hap"; done < "$out"/multiline.tsv > "$out"/multiline.txt
multiline_files="$(cat "$out"/multiline.txt | ( [ "$out" == "out" ] && cat || sed 's+out/+'"$out"'/+'))"

# let's report a few statistics
num_tot_regions="$(ls -d "$out"/*_*-* | wc -l)"
num_regions="$(echo "$multiline_files" | wc -l)"
echo "Out of $num_tot_regions regions, $num_regions has at least one haplotype with more than one variant."
echo "Of those $num_regions, here is a breakdown of the number of haplotypes each region had:"
while read hap; do grep '^H' $hap | wc -l; done < <(echo "$multiline_files") | sort | uniq -c
echo "Of those $num_regions, here is a breakdown of the number of alleles in each haplotype:"
avg_num_alleles="$(for i in $multiline_files; do grep '^V' $i | cut -f2 | sort | uniq -c | sed 's/^ *//' | cut -f1 -d' '; done | tee >(sort | uniq -c 2>&1) | awk '{ total += $1 } END { print total/NR }')"
echo "Of those $num_regions, the average number of alleles in each haplotype is $avg_num_alleles."

if [ "$mode" == "geuvadis" ]; then
  # now, let's make the variance_explained.png plot
  workflow/scripts/variance_explained_plot.py --verbosity WARNING -o "$out/variance_explained.png" -s <(echo "$multiline_files") "$out"/{locus}/happler/run/{gene}/include/merged.pgen data/geuvadis/phenos/{gene}.pheno "$out"/{locus}/happler/run/{gene}/happler.hap
elif [ "$mode" == "ukb" ]; then
  # now, let's make the variance_explained.png plot
  workflow/scripts/variance_explained_plot.py --verbosity WARNING -o "$out/variance_explained.png" -s <(echo "$multiline_files") "$out"/{locus}/happler/run/{gene}/include/merged.pgen data/ukb/phenos/{gene}.resid.pheno "$out"/{locus}/happler/run/{gene}/happler.hap
elif [ "$mode" == "aou" ]; then
  workflow/scripts/variance_explained_plot.py --verbosity WARNING -o "$out/variance_explained.png" -s <(echo "$multiline_files") "$out"/{locus}/happler/run/{gene}/include/merged.pgen data/aou/phenos/{gene}.resid.pheno "$out"/{locus}/happler/run/{gene}/happler.hap
fi
echo "Created $out/variance_explained.png" 1>&2

cd "$out"

# now, let's make the pips.tsv and exclude_pips.tsv files
echo -e "locus\tpip" > pips.tsv
for i in $(sed 's+happler.hap$+include/susie_pips.tsv+;s+out/++' multiline.txt); do grep -P '^H[0-9]\t' $i | sed 's+^+'"$(echo $i | sed 's\/include/susie_pips.tsv$\\;s+/happler/run/+:+')"':+'; done >> pips.tsv
echo "Created $out/pips.tsv" 1>&2
echo "$(awk -F $'\t' '$2 > 0.9' pips.tsv | wc -l) of the $(cat pips.tsv | wc -l) haplotypes have PIPs above 0.9"
# figure out the best PIP among only the SNPs (excluding the hap)
echo -e "locus\tpip" > exclude_pips.tsv
for i in $(sed 's+happler.hap$+exclude/susie_pips.tsv+;s+out/++' multiline.txt); do awk '(NR==1) || ($2 > max){max=$2; rec=$0} END{if (NR) print rec}' "$i" | cut -f2 | sed 's+^+'"$(echo $i | sed 's\/exclude/susie_pips.tsv$\\;s\^out/\\;s+/happler/run/+:+')"'\t+'; done >> exclude_pips.tsv
echo "Created $out/exclude_pips.tsv" 1>&2

# now, let's make the hap_pips.png file
(
  echo "a=["$(cut -f2 pips.tsv | tail -n+2 | paste -s -d,)"]"
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
echo "Created $out/hap_pips.png" 1>&2

# let's make a plot to show the haplotype PIPs vs best SuSiE PIPs when the hap is excluded
(
  echo 'a=['$(for i in $(sed 's+happler.hap$+exclude/susie_pips.tsv+;s+out/++' multiline.txt); do echo "$(awk '$1 == "H0"' "$(echo "$i" | sed 's/exclude/include/')" | cut -f2),$(awk '(NR==1) || ($2 > max){max=$2; rec=$0} END{if (NR) print rec}' "$i" | cut -f2)"; done | sed 's/^/(/;s/$/)/' | paste -s -d,)']'
  cat <<'EOF'
import numpy as np
import matplotlib.pyplot as plt
data = np.array(a)
plt.scatter(data[:,1], data[:,0]/data[:,1])
plt.axline([0, 1], [1, 1])
plt.xlabel("Best SNP PIP (when haplotype is excluded)")
plt.ylabel("Haplotype PIP / Best SNP PIP")
plt.savefig("in_vs_ex_pips.png")
EOF
) | python
echo "Created $out/in_vs_ex_pips.png" 1>&2

# now, let's check the HWE of the haplotypes
for i in $(sed 's/.hap$/.pvar/;s+out/++' multiline.txt); do plink2 --pfile ${i%.*} --hardy --out ${i%.*}-hwe &>/dev/null; done
(
  echo "a=["$(for i in $(sed 's+happler.hap$+happler-hwe.hardy+;s+out/++' multiline.txt); do cut -f10 $i | tail -n+2; done | paste -s -d,)"]"
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
echo "Created $out/hwe.png" 1>&2

# now, let's threshold by MAC and create a histogram
for i in $(sed 's/.hap$/.pvar/;s+out/++' multiline.txt); do plink2 --pfile ${i%.*} --mac 20 --make-pgen --out ${i%.*}-maf --freq &>/dev/null; done
echo "$(grep 'Error: No variants remaining' $(sed 's+out/++;s+.hap$+-maf.log+' multiline.txt) | cut -d '/' -f2,5 | wc -l) haplotypes had an MAC below 20."
(
  echo 'a=['$(cat $(sed 's+out/++;s+.hap$+-maf.afreq+' multiline.txt) | grep -Ev '^#' | cut -f 5 | paste -s -d,)']'
  cat <<'EOF'
import numpy as np
import matplotlib.pyplot as plt
data = np.array(a)
data = np.min(np.array([data, 1-data]), axis=0)
binwidth = 0.025
plt.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth))
plt.title("Haplotype MAFs")
plt.savefig("mafs.png")
EOF
) | python
echo "Created $out/mafs.png" 1>&2
echo -e "locus\tmaf" > mafs.tsv
for i in $(sed 's+out/++;s+.hap$+-maf.log+' multiline.txt); do grep -Ev '^#' $i | cut -f2,5 | sed 's+^+'"$(echo $i | sed 's\/happler-maf.afreq$\\;s\^'"$out"'/\\;s+/happler/run/+:+')"':+'; done >> mafs.tsv
echo "Created $out/mafs.tsv" 1>&2

if [ "$mode" == "geuvadis" ]; then
  # create SV LD plot
  # first, copy all of the results over
  mkdir -p sv_ld/H0
  for i in $(sed 's+out/++;s+.hap$+_svs.ld+' multiline.txt); do region="$(echo "$i" | sed 's\/happler_svs.ld$\\;s+/happler/run/+\t+')"; cp "$i" sv_ld/H0/$(echo "$region" | cut -f1):$(echo "$region" | cut -f2).ld; done
  # now, collate the results
  { echo -e 'file\tpip\tpos\tid\tld'; tail -n+2 pips.tsv | sort -gr -k2,2 | { while read -r line; do file="sv_ld/$(echo "$line" | cut -f1 | cut -f3 -d:)/$(echo "$line" | cut -f1 | cut -f-2 -d:).ld"; echo -en "$file"$'\t'; echo -en "$(echo "$line" | cut -f2)"$'\t'; awk -F $'\t' -v 'OFS=\t' '{print $2, $3, sqrt($4*$4);}' "$file" | sort -gr -k3,3 | head -n1; done } | sed 's/^.*sv_ld\///'; } > pips_sv_ld.tsv
  echo "Created $out/pips_sv_ld.tsv" 1>&2
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
  echo "Created $out/pips_best_sv_ld.png" 1>&2

  cd -
  # create STR LD plot
  # first, let's get a list of the STRs in the regions with haplotypes
  (
    echo -ne "region\tgene\thap.id\tpip\t"
    head -n1 data/geuvadis/mlamkin/Geuvadis_varlevel_corrected_significant_variants.with-end.tsv && \
    ~/miniconda3/envs/htslib/bin/bedtools intersect -a <(
      echo -e "chrom\tstart\tend\tgene\thap.id\tpip" && tail -n+2 "$out"/pips.tsv | sed 's+_+\t+;s+-+\t+;s+:+\t+g' | sort -k1,1V -k2,2n
    ) -b data/geuvadis/mlamkin/Geuvadis_varlevel_corrected_significant_variants.with-end.tsv -wa -wb -loj | \
    sed 's+\t+_+;s+\t+-+'
  ) | awk -F'\t' '$2 == $8' | cut -f8 --complement > "$out"/STR_assocations.tsv
  echo "Created $out/STR_assocations.tsv" 1>&2
  # now, let's compute LD for each region
  echo -e "hap\tpip\tpos\tid\tld\talleles" > "$out"/pips_str_ld.tsv
  while IFS= read -r line; do
    str_id="$(echo "$line" | cut -f5,6 --output-delimiter ':')"
    echo -ne "$(echo "$line" | cut -f1-3 --output-delimiter ':')\t$(echo "$line" | cut -f4)\t$(echo "$line" | cut -f6)\t$str_id\t"
    echo -ne "$(workflow/scripts/compute_pgen_ld.py --verbosity WARNING --target-is-repeat --hap-id "$str_id" -o /dev/stdout "$out/$(echo "$line" | cut -f1)"/happler/run/"$(echo "$line" | cut -f2)"/happler.pgen data/geuvadis/mlamkin/all_Geuvadis_STRs.pgen | tail -n+2 | cut -f4)"
    echo -e "\t$(grep -P '\t'"$str_id"'\t' data/geuvadis/mlamkin/all_Geuvadis_STRs.pvar | cut -f 4,5 --output-delimiter ,)"
  done < <(tail -n+2 "$out"/STR_assocations.tsv) >> "$out"/pips_str_ld.tsv
  echo "Created $out/pips_str_ld.tsv" 1>&2
  (
    head -n1 "$out/pips_str_ld.tsv"
    tail -n+2 "$out/pips_str_ld.tsv" \
      | awk -F'\t' -v OFS='\t' '{$5 = ($5 < 0) ? -$5 : $5; print}' \
      | sort -t$'\t' -k1,1 -k5,5nr \
      | awk -F'\t' -v OFS='\t' '!seen[$1]++ { print }'
  ) > "$out/pips_best_str_ld.tsv"
  echo "Created $out/pips_best_str_ld.tsv" 1>&2
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
  echo "Created $out/pips_best_str_ld.png" 1>&2
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
  echo "Created $out/ld_str_vs_sv.png"; 1>&2
fi

# let's make a plot to show runtime and memory usage
echo -e "locus\tnum_vars\ttime_s\tmem_mb" > bench.tsv
# To get just the multiline ones, use '$(sed 's+out/++;s+happler.hap$+bench/run+' multiline.txt)' instead of */happler/run/*/bench/run
for i in */happler/run/*/bench/run; do
  echo -e "$(echo $i | sed 's+/happler/run/+:+;s+/bench/run++')\t$(wc -l "$(echo "$i" | sed 's+happler/run+genotypes+;s+bench/run+snps.pvar+')" | cut -f1 -d' ')\t$(cut -f1,3 "$i" | tail -n1)"
done >> bench.tsv
echo "Created $out/bench.tsv" 1>&2
(
  echo "a=["$(tail -n+2 bench.tsv | cut -f2- --output-delimiter , | sed 's+^+(+;s+$+)+' | paste -s -d,)"]"
  cat <<'EOF'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = np.array(a)
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(6.5,3))
axes[0].scatter(data[:, 0], data[:, 1]/60)
axes[0].set_xlabel("Number of variants in locus")
axes[0].set_ylabel("Happler Runtime (mins)")
axes[1].scatter(data[:, 0], data[:, 2]/1000)
axes[1].set_xlabel("Number of variants in locus")
axes[1].set_ylabel("Happler Max Memory (GB)")
plt.tight_layout()
plt.savefig("bench.png")
EOF
) | python
echo "Created $out/bench.png" 1>&2

# merge all of the .tsv files together
(
  echo -e "locus\thap_pip\tbest_snp_pip"
  join -t $'\t' -j1 <(tail -n+2 pips.tsv | sed 's/:/\t/g' | sed 's/\t/:/' | sort -k1,1) <(tail -n+2 exclude_pips.tsv | sort -k1,1) | sed 's/\t/:/' | sort -k1,1
) > merged_pips.tsv

join -t $'\t' -j1 --header merged_pips.tsv <(
  head -n1 variance_explained.tsv
  tail -n+2 variance_explained.tsv | sort -k1,1
) | sort -k1,1 | (
  if [ "$mode" == "geuvadis" ]; then
    # also merge with the last four columns of pips_best_str_ld.tsv and the last three columns of pips_sv_ld.tsv
    join -t $'\t' -j1 - <(
      join -t $'\t' -j1 --header <(
        echo -ne "locus\t"; head -n1 pips_sv_ld.tsv | cut -f3- | tr $'\t' $'\n' | sed 's/^/sv_/' | paste -s
        tail -n+2 pips_sv_ld.tsv | cut -f1,3- | sed 's+/+~+;s+.ld\t+~+' | awk -F '~' -v 'OFS=\t' '{print $2":"$1,$3;}' | sort -k1,1
      ) <(
        echo -ne "locus\t"; head -n1 pips_best_str_ld.tsv | cut -f3- | tr $'\t' $'\n' | sed 's/^/str_/' | paste -s
        tail -n+2 pips_best_str_ld.tsv | cut -f1,3- | sort -k1,1
      ) | sort -k1,1
    )
  else
    cat
  fi
) | sort -k1,1g > merged.tsv
echo "Created $out/merged.tsv" 1>&2

if [ "$mode" == "aou" ]; then
  zip multiline.zip variance_explained.png hap_pips.png in_vs_ex_pips.png hwe.png mafs.png bench.png merged.tsv
  echo "Created $out/multiline.zip" 1>&2
fi
