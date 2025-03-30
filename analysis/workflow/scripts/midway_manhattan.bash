#!/usr/bin/env bash
# Creates a manhattan plot of SNP p-values midway through happler's tree branching process
# arg1: PGEN file from which to obtain SNP genotypes
# arg2: pheno file from which to obtain the phenotypes
# arg3: happler hap file from which to obtain the haplotypes
# arg4: output prefix
# arg5: ID of the target haplotype in the hap file
# arg6: ID of a "child" SNP in the target haplotype. The parent will include all SNPs of the haplotype up to this one.
# arg7: MAF threshold for filtering SNPs
# arg8: 0 or 1 indicating whether to include just the parent (0) or both the parent and child SNPs (1) in the regression as covariates, or 2 if the regular p-values should be converted to t-test p-values. 1 requires that the child SNP be the second node from the root of the tree. Additionally, the value 3 indicates that t-test p-values should be computed with a covariance correction and the value 4 indicates that delta BIC values should be generated in place of p-values. (optional - defaults to 0)
# arg9: 0 or 1 indicating whether to compute p-values for just the target SNP or all SNPs.  (optional - defaults to 0)
# ex: workflow/scripts/midway_manhattan.bash out/19_55363180-55833573/genotypes/snps.pgen 19_55363180-55833573.indep.pheno 19_55363180-55833573.indep.hap 19_55363180-55833573.indep H0 rs61734259 0.005 1

pgen_file="$1"
pheno_file="$2"
hap_file="$3"
out_prefix="$4"
hap_id="$5"
snp_id="$6"
maf="$7"
condition="${8:-0}"
just_target="${9:-0}"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

############################################## CHECKS ##############################################

# first, check that the child SNP is in the hap file
grep -m1 -P '^V\t'"$hap_id"'\t.*\t'"$snp_id"'\t' "$hap_file" >/dev/null || {
    echo "Failed to find $hap_id:$snp_id in $hap_file" 1>&2
    exit
}

# now, retrieve the contents of the hap file that denote the haplotype up until (but not including!) the child SNP
hap="$(grep -P '\t'"$hap_id"'\t' "$hap_file" | sed '/'"$snp_id"'/Q')"

# next, check that the haplotype has at least one SNP preceding the child SNP
[ "$(echo "$hap" | wc -l)" -gt 1 ] || {
    echo "Failed to find at least one parent SNP for $hap_id:$snp_id in $hap_file" 1>&2
    exit
}

# if we should include the parent as a covariate, then there should ONLY be one SNP preceding the child SNP
if [ "$condition" -eq 1 ]; then
    [ "$(echo "$hap" | wc -l)" -eq 2 ] || {
        echo "There can only be one SNP before $snp_id for $hap_id in $hap_file" 1>&2
        exit
    }
fi
parent_snp_id="$(echo "$hap" | tail -n1 | cut -f5)"

# finally, check that the child SNP appears in the PVAR file
grep -m1 -P "\t"$snp_id"\t" "${pgen_file%.pgen}.pvar" >/dev/null || {
    echo "Failed to find $hap_id:$snp_id in ${pgen_file%.pgen}.pvar" 1>&2
    exit
}

# now, determine whether the child SNP appears in the haplotype by its REF (0) or ALT (1) allele
allele=$(grep "$(grep -m1 -P '^V\t'"$hap_id"'\t.*\t'"$snp_id"'\t' "$hap_file" | cut -f5,6)" "${pgen_file%.pgen}.pvar" | wc -l)
allele=$(expr 1 - $allele)

# use --mac 1 if maf is 0
[ "$maf" == "0" ] && maf="0.0000000005"

# set up just_target to contain the proper parameter
if [ "$just_target" == "1" ]; then
    just_target="--var-id $snp_id"
else
    just_target=""
fi

############################################## MAIN PROGRAM ########################################

# step 1: use happler transform to retrieve the genotypes of the parent haplotype extended by each SNP in the PGEN file
happler transform \
--maf $maf \
-a $allele \
--verbosity DEBUG \
-o "$out_prefix".pgen \
-S <(cut -f1 "$pheno_file" | tail -n+2) \
$just_target \
"$pgen_file" <(
    grep -E '^#' "$hap_file"
    echo "$hap"
)

# check that the child SNP is in the transformed output
grep -m1 -P "\t"$snp_id"\t" "$out_prefix".pvar >/dev/null || {
    echo -e "Failed to find $hap_id:$snp_id in $out_prefix.pvar\nCheck the haplotype's MAF" 1>&2
    exit
}

# step 2: use plink2 to obtain p-values for a manhattan plot
if [ "$condition" -eq 1 ]; then
    # output the parent and child nodes to a covariate file
    plink2 \
    --maf $maf \
    --nonfounders \
    --out "$out_prefix" \
    --export A ref-first \
    --pfile "${pgen_file%.pgen}" \
    --snps $parent_snp_id $snp_id
    { echo -n "#"; cut -f 2,7- "$out_prefix".raw; } > "$out_prefix".covar
    # incorporate the parent and child nodes as covariates in the model
    plink2 \
    --maf $maf \
    --nonfounders \
    --vif 10000000000 \
    --max-corr 0.9999999999 \
    --glm no-x-sex hide-covar \
    --out "$out_prefix" \
    --pfile "$out_prefix" \
    --variance-standardize \
    --no-input-missing-phenotype \
    --pheno iid-only "$pheno_file" \
    --covar 'iid-only' "$out_prefix".covar

    last_arg="-a 0.05"
else
    plink2 \
    --maf $maf \
    --nonfounders \
    --out "$out_prefix" \
    --pfile "$out_prefix" \
    --variance-standardize \
    --no-input-missing-phenotype \
    --pheno iid-only "$pheno_file" \
    --glm no-x-sex allow-no-covars

    last_arg="-a $(grep -P '^V\t.*\t'"$parent_snp_id"'\t' "$hap_file" | cut -f7)"
fi

linear_file="$(ls -1 "$out_prefix".*.glm.linear | grep -v 'out\.parent\..*\.glm\.linear' | tail -n1)"

if [ "$condition" -eq 2 ] || [ "$condition" -eq 3 ] || [ "$condition" -eq 4 ]; then
    # first, get the parent haplotype as a PGEN file
    haptools transform \
    --verbosity DEBUG \
    -o "$out_prefix".parent.pgen \
    -S <(cut -f1 "$pheno_file" | tail -n+2) \
    "$pgen_file" <(
        grep -E '^#' "$hap_file"
        echo "$hap"
    )
    # get summary statistics for the parent haplotype
    plink2 \
    --maf $maf \
    --nonfounders \
    --variance-standardize \
    --out "$out_prefix".parent \
    --pfile "$out_prefix".parent \
    --no-input-missing-phenotype \
    --pheno iid-only "$pheno_file" \
    --glm no-x-sex allow-no-covars
    # set up the --bic or --covariance flag, if needed
    extra_flag=""
    if [ "$condition" -eq 3 ]; then
        extra_flag="--covariance"
    elif [ "$condition" -eq 4 ]; then
        extra_flag="--bic"
    fi
    # now, convert the plink2 --glm results into t-test p-values
    "$SCRIPT_DIR"/linear2ttest.py \
    --bic \
    -i "$hap_id" \
    --verbosity DEBUG \
    -o "$out_prefix".ttest.linear \
    "$linear_file" "$out_prefix".parent.*.glm.linear \
    "$out_prefix".pgen "$out_prefix".parent.pgen \
    "$pheno_file"
    # rename the linear file
    linear_file="$out_prefix".ttest.linear
    last_arg="-a 0.05"
fi

# save symlink to linear file for future use
ln -sfnr "$linear_file" "$out_prefix".linear

if [ "$just_target" == "" ]; then
    # step 3: use manhattan.py to actually create the manhattan plot
    "$SCRIPT_DIR"/manhattan.py \
    --bic \
    -i "$snp_id" \
    -o "$out_prefix".png \
    -l "$linear_file" \
    $last_arg
fi
