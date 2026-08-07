#!/usr/bin/env bash

# arg1: the output path from AoU (ex: out)
# arg2: MAFs (ex: "0.001 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.065 0.085 0.1")
# ex: workflow/scripts/maf-figures.bash out "0.001 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.065 0.085 0.1"

out="$1"
maf_str="$2"
read -r -a mafs <<< "$maf_str"

cd "$out/pip_maf"

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

for f in glob.glob("[0-9]*.tsv"):
  f_tab = np.genfromtxt(fname=f, dtype=None, delimiter="\t", names=True, encoding='utf-8')
  if not len(f_tab) or f_tab['snp_maf'][0] > max(mafs):
    continue
  for hap_id in np.unique(f_tab['hap_id']):
    diffs_hap = {maf: np.nan for maf in mafs}
    f_tab_hap = f_tab[f_tab['hap_id'] == hap_id]
    diffs_hap.update(dict(zip(f_tab_hap['maf_thresh'], f_tab_hap['hap_r2']-f_tab_hap['snp_r2'])))

    diffs.append(diffs_hap)
    labels.append(f.removesuffix(".tsv").replace(" ","\t").replace('\uf03a', ':')+"\t"+hap_id)

data = np.array([list(diff_hap.values()) for diff_hap in diffs])

labels_arr = np.array(labels)
combined = np.empty((data.shape[0], data.shape[1] + 1), dtype=object)
combined[:, 0] = labels_arr
combined[:, 1:] = data

header_str = 'locus\tcausal_snp\thap_id\t' + '\t'.join(map(str, mafs))
formats = ['%s'] + ['%f'] * data.shape[1]

np.savetxt('hap_snp_r2.tab', combined, delimiter='\t', header=header_str, comments='', fmt=formats)

# drop any rows where <= 1 value is not np.nan
data = data[np.sum(~np.isnan(data), axis=1) > 1]
# combined = combined[np.sum(~np.isnan(data), axis=1) > 1]

plt.figure(figsize=(10, 6.5))
for row in data:
    plt.plot(mafs, row, marker='o', alpha=0.3)
plt.xlabel('MAF Threshold')
plt.ylabel('Difference (hap_r2 - snp_r2)')
plt.savefig('hap_snp_r2_diff.full.pdf')

# drop any rows with negative values
data = data[~np.any(data < 0, axis=1)]

plt.clf()
plt.figure(figsize=(10, 6.5))
for row in data:
    plt.plot(mafs, row, marker='o', alpha=0.3)
plt.xlabel('MAF Threshold')
plt.ylabel('Difference (hap_r2 - snp_r2)')
plt.savefig('hap_snp_r2_diff.pdf')
EOF
) | python
echo "Created $out/ld_maf/hap_snp_r2.tab" 1>&2
echo "Created $out/ld_maf/hap_snp_r2_diff.pdf" 1>&2
echo "Created $out/ld_maf/hap_snp_r2_diff.full.pdf" 1>&2
