from pathlib import Path


out = config["out"] + "/genotypes"
logs = out + "/logs"
bench = out + "/bench"


locus_chr = config["locus"].split(":")[0]
locus_start = config["locus"].split(":")[1].split('-')[0]
locus_end = config["locus"].split(":")[1].split('-')[1]


def check_config(value, default=False, place=config, as_set=False):
    """return true if config value exists and is true"""
    value = place[value] if (value in place and place[value]) else default
    return (set(value) if isinstance(value, list) else {value}) if as_set else value


rule plink2vcf:
    """ convert a PLINK file to VCF """
    input:
        pgen=lambda wildcards: expand(config["snp_panel"], chr=locus_chr)[0],
    params:
        pfile=lambda wildcards, input: str(Path(input.pgen).with_suffix('')),
        out=lambda wildcards, output: str(Path(output.bcf).with_suffix('')),
        start=locus_start,
        end=locus_end,
        chrom=locus_chr,
    output:
        vcf=temp(out + "/unphased.vcf.gz"),
        bcf=out + "/unphased.bcf",
        bcf_idx=out + "/unphased.bcf.csi",
        log=temp(out + "/unphased.log"),
    resources:
        runtime="0:20:00"
    log:
        logs + "/plink2vcf"
    benchmark:
        bench + "/plink2vcf"
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --pfile {params.pfile} --out {params.out} --from-bp {params.start} "
        "--to-bp {params.end} --chr {params.chrom} --snps-only 'just-acgt' "
        "--export vcf bgz id-paste=iid &>{log} && "
        "bcftools view -Ob -o {output.bcf} {output.vcf} &>>{log} && "
        "tabix -p bcf {output.bcf} &>>{log}"


def phase_gt_input(wildcards):
    input_files = { 'map': config['phase_map'].format(chr=locus_chr) }
    if config["snp_panel"].endswith('.pgen'):
        input_files['unphased'] = rules.plink2vcf.output.bcf
    else:
        input_files['unphased'] = config["snp_panel"]
    input_files['unphased_idx'] = input_files['unphased'] + ".csi"
    return input_files


rule phase_gt:
    """ phase an unphased set of genotypes """
    input: unpack(phase_gt_input)
    params:
        locus=config["locus"],
    output:
        phased=out + "/phased.bcf",
        phased_idx=out + "/phased.bcf.csi",
    resources:
        runtime="4:00:00"
    threads: 8
    log:
        logs + "/phase_gt"
    benchmark:
        bench + "/phase_gt"
    conda:
        "../envs/shapeit.yml"
    shell:
        "shapeit4 --input {input.unphased} --map {input.map} --region {params.locus} "
        "--output {output.phased} --thread {threads} --log {log} &>>{log} && "
        "tabix -p bcf {output.phased} &>>{log}"


rule keep_samps:
    """ create a list of the samples that we should keep """
    input:
        snp_vcf=lambda wildcards: rules.phase_gt.output.phased if check_config('phase_map')
            else config["snp_panel"],
        snp_vcf_idx=lambda wildcards: rules.phase_gt.output.phased_idx if check_config('phase_map')
            else config["snp_panel"] + ".tbi",
        str_vcf=expand(config["str_panel"], chr=locus_chr),
        str_vcf_idx=expand(config["str_panel"], chr=locus_chr),
        samp=lambda wildcards: config["exclude_samples"],
    output:
        samples=out+"/samples.tsv"
    log:
        logs+"/keep_samps"
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/keep_samps.bash {input.snp_vcf} {input.str_vcf} {input.samp}"
        " > {output.samples} 2>{log}"


rule vcf2plink:
    """ convert a VCF into a PGEN file """
    input:
        vcf=lambda wildcards: rules.phase_gt.output.phased if check_config('phase_map')
            else config["snp_panel"],
        vcf_idx=lambda wildcards: rules.phase_gt.output.phased_idx if check_config('phase_map')
            else config["snp_panel"] + ".tbi",
        samples=rules.keep_samps.output.samples if check_config('exclude_samples') else [],
    params:
        maf=config["min_maf"],
        prefix=lambda wildcards, output: Path(output.pgen).with_suffix(""),
    output:
        pgen=out+"/snps.pgen",
        pvar=out+"/snps.pvar",
        psam=out+"/snps.psam",
        log=temp(out+"/snps.log"),
    resources:
        runtime="0:03:00"
    log:
        logs + "/vcf2plink",
    benchmark:
        bench + "/vcf2plink",
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --bcf {input.vcf} --maf {params.maf} --make-pgen "
        "--keep {input.samples} " if check_config('exclude_samples') else ""
        " --out {params.prefix} &>{log}"


rule subset_str:
    """ subset samples from a STR VCF """
    input:
        vcf=lambda wildcards: expand(config["str_panel"], chr=locus_chr)[0],
        vcf_idx=lambda wildcards: expand(config["str_panel"] + ".tbi", chr=locus_chr),
        samples=rules.keep_samps.output.samples,
    output:
        vcf=out+"/strs.bcf",
        vcf_idx=out+"/strs.bcf.csi",
    resources:
        runtime="0:30:00"
    log:
        logs + "/subset_str",
    benchmark:
        bench + "/subset_str",
    conda:
        "../envs/default.yml"
    shell:
        "bcftools view -S {input.samples} --write-index "
        "-Ob -o {output.vcf} {input.vcf}"
