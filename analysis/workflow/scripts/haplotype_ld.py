#!/usr/bin/env python
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


sns.set(rc={'figure.figsize':(11.7,8.27)})

phens = pd.read_csv("../phens.tsv.gz", sep="\t", index_col=0, header=0)
gens = pd.read_csv("genotypes.tsv", sep="\t", index_col=1)
gt = gens.iloc[:,3:]
alleles = gens['alleles']

SNPs = gt.iloc[~gt.index.str.startswith('STR')]
STR = gt.iloc[gt.index.str.startswith('STR')]

samples = phens.index.astype(str)
num_samp = len(samples)

lens = {str(i): len(al) for i, al in enumerate(gens.loc[STR.index[0],'alleles'].split(','))}
str_gt = STR.iloc[0].str.split('|', expand=True).applymap(lens.__getitem__)
str_gt = str_gt.loc[samples]
str_gts = pd.concat((str_gt[0], str_gt[1]), axis=0, ignore_index=True)
split_chroms = pd.DataFrame({STR.index[0]: str_gts})

fig, axes = plt.subplots(1, len(SNPs.index))
for idx, snp in enumerate(SNPs.index):
    snp_gt = gt.loc[snp].str.split('|', expand=True)
    snp_gt = snp_gt.loc[samples]
    snp_gts = pd.concat((snp_gt[0], snp_gt[1]), axis=0, ignore_index=True).astype(int)
    split_chroms[snp] = snp_gts
    sns.violinplot(ax=axes[idx], x=snp, y=STR.index[0], data=split_chroms, order=[0,1])

plt.savefig('snps.png')
plt.clf()

fig, axes = plt.subplots(1, 2)
hap_alleles = {snp: 1 for snp in SNPs.index}
hap = pd.DataFrame([(split_chroms[snp] == hap_alleles[snp]).astype(bool) for snp in SNPs.index]).T
hap['haplotype'] = hap.all(axis=1).astype(int)
hap[STR.index[0]] = split_chroms[STR.index[0]]
sns.violinplot(ax=axes[0], x='haplotype', y=STR.index[0], data=hap, order=[0,1])

hap_sum = pd.DataFrame({
    'haplotype': (hap['haplotype'].iloc[:num_samp].to_numpy() + hap['haplotype'].iloc[num_samp:].to_numpy()),
    STR.index[0]: (hap[STR.index[0]].iloc[:num_samp].to_numpy() + hap[STR.index[0]].iloc[num_samp:].to_numpy()),
})
sns.regplot(ax=axes[1], x='haplotype', y=STR.index[0], data=hap_sum)
plt.savefig('haplotype.png')
