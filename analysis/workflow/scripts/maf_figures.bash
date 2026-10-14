#!/usr/bin/env bash

# arg1: the output path from AoU (ex: out)
# arg2: MAFs (ex: "0.001 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.065 0.085 0.1")
# ex: workflow/scripts/maf-figures.bash out "0.001 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.065 0.085 0.1"

out="$1"
maf_str="$2"
read -r -a mafs <<< "$maf_str"

cd "$out/pip_maf"
(
  echo "mafs=["$(echo "${mafs[@]}" | sed 's/ /,/g')"]"
  cat <<'EOF'
import glob
import matplotlib
import numpy as np
matplotlib.use('Agg')
from pathlib import Path
import matplotlib.pyplot as plt

diffs = []
labels = [] # List to track the file names

for f in glob.glob("[0-9]*.tsv"):
  f_tab = np.genfromtxt(fname=f, dtype=None, delimiter="\t", names=True, encoding='utf-8')
  for hap_id in np.unique(f_tab['hap_id']):
    diffs_hap = {maf: np.nan for maf in mafs}
    f_tab_hap = f_tab[f_tab['hap_id'] == hap_id]
    diffs_hap.update(dict(zip(f_tab_hap['maf_thresh'], f_tab_hap['hap_pip']-f_tab_hap['snp_pip'])))

    diffs.append(diffs_hap)
    labels.append(f.removesuffix(".tsv").replace(" ","\t").replace('\uf03a', ':')+"\t"+hap_id)

data = np.array([list(diff_hap.values()) for diff_hap in diffs])

labels_arr = np.array(labels)
combined = np.empty((data.shape[0], data.shape[1] + 1), dtype=object)
combined[:, 0] = labels_arr
combined[:, 1:] = data

header_str = 'locus\thap_id\t' + '\t'.join(map(str, mafs))
formats = ['%s'] + ['%f'] * data.shape[1]

np.savetxt('hap_snp_pip.tab', combined, delimiter='\t', header=header_str, comments='', fmt=formats)

# drop any rows where <= 1 value is not np.nan
data = data[np.sum(~np.isnan(data), axis=1) > 1]
# combined = combined[np.sum(~np.isnan(data), axis=1) > 1]

plt.figure(figsize=(10, 6.5))
for row in data:
    plt.plot(mafs, row, marker='o', alpha=0.3)
plt.xlabel('MAF Threshold')
plt.ylabel('Difference (hap_pip - snp_pip)')
plt.savefig('hap_snp_pip_diff.pdf')
EOF
) | python
echo "Created $out/pip_maf/hap_snp_pip.tab" 1>&2
echo "Created $out/pip_maf/hap_snp_pip_diff.pdf" 1>&2
cd -
cd "$out/ld_maf"

for i in *.tsv; do head -n2 "$i" | tail -n1 | cut -f5; done > rare_snp_mafs.txt
echo "Created $out/ld_maf/rare_snp_mafs.txt" 1>&2

# now, let's make the rare_snp_mafs.pdf file
(
  echo "a=["$(cat rare_snp_mafs.txt | paste -s -d,)"]"
  cat <<'EOF'
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = np.array(a)
data = -np.log10(data[data < 0.1])
binwidth = 0.01
plt.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth))
plt.title("Causal variant MAFs")
plt.tight_layout()
plt.savefig("rare_snp_mafs.pdf", bbox_inches="tight")
EOF
) | python
echo "Created $out/ld_maf/rare_snp_mafs.pdf" 1>&2

(
  echo "mafs=["$(echo "${mafs[@]}" | sed 's/ /,/g')"]"
  cat <<'EOF'
import glob
import matplotlib
import numpy as np
matplotlib.use('Agg')
from pathlib import Path
import matplotlib.pyplot as plt

diffs = []
labels = [] # List to track the file names
ignore = []
low_delta_bic = ["14_67259065-69458750", "2_158835987-161011584"]
high_pip = []

pips_mac20 = np.genfromtxt(fname="../loci_best_snp_pip_mac20.tsv", dtype=None, delimiter="\t", names=["locus","snp_id","pip"], encoding='utf-8')
pips_mac20 = ["\t".join(i) for i in list(pips_mac20[pips_mac20["pip"] > 0.9][["locus", "snp_id"]])]

for f in glob.glob("[0-9]*.tsv"):
  f_tab = np.genfromtxt(fname=f, dtype=None, delimiter="\t", names=True, encoding='utf-8')
  label = f.removesuffix(".tsv").replace(" ","\t").replace('\uf03a', ':')
  if not len(f_tab) or f_tab['snp_maf'][0] > max(mafs) or label.split("\t")[0] in low_delta_bic:
    continue
  for hap_id in np.unique(f_tab['hap_id']):
    diffs_hap = {maf: np.nan for maf in mafs}
    f_tab_hap = f_tab[f_tab['hap_id'] == hap_id]
    diffs_hap.update(dict(zip(f_tab_hap['maf_thresh'], f_tab_hap['hap_r2']-f_tab_hap['snp_r2'])))

    if len(f_tab_hap['snp_maf']) > 1:
      next_lowest_maf = np.unique(f_tab_hap['snp_maf'])[1]
      snp_r2_next_lowest_maf = np.unique(f_tab_hap[f_tab_hap['snp_maf'] == next_lowest_maf]['snp_r2'])[0]
    else:
      snp_r2_next_lowest_maf = float("inf")

    label_hap = label+"\t"+hap_id

    diffs.append(diffs_hap)
    labels.append(label_hap)
    if snp_r2_next_lowest_maf < 0.05:
      ignore.append(label_hap)
    if label in pips_mac20:
      high_pip.append(label_hap)

data = np.array([list(diff_hap.values()) for diff_hap in diffs])

labels_arr = np.array(labels)
combined = np.empty((data.shape[0], data.shape[1] + 1), dtype=object)
combined[:, 0] = labels_arr
combined[:, 1:] = data

header_str = 'locus\tcausal_snp\thap_id\t' + '\t'.join(map(str, mafs))
formats = ['%s'] + ['%f'] * data.shape[1]

np.savetxt('hap_snp_r2.tab', combined, delimiter='\t', header=header_str, comments='', fmt=formats)

unique_loci = np.unique([label.split('\t')[0] for label in combined[:, 0]])
cmap = plt.get_cmap('tab20')
# Using modulo 20 ensures we don't run out of colors if there are >20 loci
locus_colors = {locus: cmap(i % 20) for i, locus in enumerate(unique_loci)}

# drop any rows where <= 1 value is not np.nan
combined = combined[np.sum(~np.isnan(data), axis=1) > 1]
data = data[np.sum(~np.isnan(data), axis=1) > 1]

plt.figure(figsize=(5.4, 3.5))
for label, row in zip(combined[:, 0], data):
    locus = label.split('\t')[0]
    plt.plot(mafs, row, marker='o', alpha=0.3, color=locus_colors[locus])
plt.xlabel('MAF Threshold')
plt.ylabel('Difference (hap_r2 - snp_r2)')
plt.savefig('hap_snp_r2_diff.full.full.pdf')
plt.close()

# drop any rows where the next best SNP wasn't even close to tagging the causal variant
data = data[~np.isin(combined[:, 0], ignore)]
combined = combined[~np.isin(combined[:, 0], ignore)]

plt.figure(figsize=(5.4, 3.5))
for label, row in zip(combined[:, 0], data):
    locus = label.split('\t')[0]
    plt.plot(mafs, row, marker='o', alpha=0.3, color=locus_colors[locus])
plt.xlabel('MAF Threshold')
plt.ylabel('Difference (hap_r2 - snp_r2)')
plt.savefig('hap_snp_r2_diff.full.pdf')
plt.close()

# drop any rows where the PIP was below 0.9 at an MAC of 20
data = data[np.isin(combined[:, 0], high_pip)]
combined = combined[np.isin(combined[:, 0], high_pip)]

plt.figure(figsize=(5.4, 3.5))
for label, row in zip(combined[:, 0], data):
    locus = label.split('\t')[0]
    plt.plot(mafs, row, marker='o', alpha=0.3)
plt.xlabel('MAF Threshold')
plt.ylabel('Difference (hap_r2 - snp_r2)')
plt.savefig('hap_snp_r2_diff.pdf')
plt.close()
EOF
) | python
echo "Created $out/ld_maf/hap_snp_r2.tab" 1>&2
echo "Created $out/ld_maf/hap_snp_r2_diff.pdf" 1>&2
echo "Created $out/ld_maf/hap_snp_r2_diff.full.pdf" 1>&2
echo "Created $out/ld_maf/hap_snp_r2_diff.full.full.pdf" 1>&2
