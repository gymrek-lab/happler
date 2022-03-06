#!/usr/bin/env bash

# arg1: the locus from which to obtain the genotypes
# arg2: VCF file from which to obtain SNP genotypes
# arg3: VCF file from which to obtain STR genotypes
# arg4: a list of samples to keep
# arg5: temporary file to write SNP VCF
# arg6: temporary file to write STR VCF
# ex: workflow/scripts/matrix.bash '1:98001984-99001984' out/1_98001984-99001984/phased.bcf strs/chr1.vcf.gz out/1_98001984-99001984/samples.tsv.gz | head


loc="$1"
snp_vcf="$2"
str_vcf="$3"
samps_file="$4"
snp_temp_file="$5"
str_temp_file="$6"


bcftools view -Oz -S <(
    zcat "$samps_file"
) -o "$temp_file" "$snp_vcf" "$loc"

bcftools view -Ou -S <(
    zcat "$samps_file" | awk '{print $1 "_" $1;}'
) "$str_vcf" "$loc" | \
bcftools reheader -o "$str_temp_file" -s <(
    zcat "$samps_file"
) -

echo -ne "POS\tID\tMAF\talleles\t"
paste -s -d $'\t' <(zcat "$samps_file")
bcftools concat -Ou "$snp_temp_file" "$str_temp_file" | {
    bcftools +fill-tags - -- -t MAF | \
    bcftools query -f '%POS\t%ID\t%INFO/MAF\t%REF,%ALT\t[%GT\t]\n' | \
    sed 's/\t$//'
}
