from pathlib import Path


out = config["out"] + "/genotypes"
logs = out + "/logs"
bench = out + "/bench"


def check_config(value, default=False, place=config, as_set=False):
    """return true if config value exists and is true"""
    value = place[value] if (value in place and place[value]) else default
    return (set(value) if isinstance(value, list) else {value}) if as_set else value


def parse_locus(locus):
    """parse locus into chrom, start, end"""
    chrom = locus.split("_")[0]
    end = locus.split("-")[1]
    start = locus.split("_")[1].split("-")[0]
    return chrom, start, end


rule plink2vcf:
    """ convert a PLINK file to VCF """
    input:
        pgen=lambda wildcards: expand(config["snp_panel"], chr=parse_locus(wildcards.locus)[0])[0],
    params:
        pfile=lambda wildcards, input: str(Path(input.pgen).with_suffix('')),
        out=lambda wildcards, output: str(Path(output.bcf).with_suffix('')),
        start=lambda wildcards: wildcards.locus.split("_")[1].split("-")[0],
        end=lambda wildcards: wildcards.locus.split("-")[1],
        chrom=lambda wilcards: wildcards.locus.split("_")[0],
    output:
        vcf=temp(out + "/unphased.vcf.gz"),
        bcf=out + "/unphased.bcf",
        bcf_idx=out + "/unphased.bcf.csi",
        log=temp(out + "/unphased.log"),
    resources:
        runtime=20,
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
    input_files = { 'map': config['phase_map'].format(chr=parse_locus(wildcards.locus)[0]) }
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
        runtime=60*4,
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
        str_vcf=lambda wildcards: expand(config["str_panel"], chr=parse_locus(wildcards.locus)[0]),
        str_vcf_idx=lambda wildcards: expand(config["str_panel"], chr=parse_locus(wildcards.locus)[0]),
        samp=lambda wildcards: config["exclude_samples"],
    output:
        samples=out+"/samples.tsv"
    resources:
        runtime=3,
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
        runtime=3,
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
        vcf=lambda wildcards: expand(config["str_panel"], chr=parse_locus(wildcards.locus)[0])[0],
        vcf_idx=lambda wildcards: expand(config["str_panel"] + ".tbi", chr=parse_locus(wildcards.locus)[0]),
        samples=rules.keep_samps.output.samples,
    output:
        vcf=out+"/strs.bcf",
        vcf_idx=out+"/strs.bcf.csi",
    resources:
        runtime=30,
    log:
        logs + "/subset_str",
    benchmark:
        bench + "/subset_str",
    conda:
        "../envs/default.yml"
    shell:
        "bcftools view -S {input.samples} --write-index "
        "-Ob -o {output.vcf} {input.vcf}"


def subset_input():
    if check_config("phase_map") or check_config("exclude_samples") or not config["snp_panel"].endswith(".pgen"):
        return {
            "pgen": rules.vcf2plink.output.pgen,
            "pvar": rules.vcf2plink.output.pvar,
            "psam": rules.vcf2plink.output.psam,
        }
    else:
        return {
            "pgen": config["snp_panel"],
            "pvar": Path(config["snp_panel"]).with_suffix(".pvar"),
            "psam": Path(config["snp_panel"]).with_suffix(".psam"),
        }


rule subset:
    """subset the simulation dataset if needed"""
    input:
        **subset_input()
    params:
        prefix=lambda wildcards, input: Path(input.pgen).with_suffix(""),
        out=lambda wildcards, output: Path(output.pgen).with_suffix(""),
        sampsize=lambda wildcards: wildcards.sampsize,
        region=lambda wildcards: locus.replace("_", ":"),
    output:
        pgen=out+"/subset/{sampsize}.pgen",
        pvar=out+"/subset/{sampsize}.pvar",
        psam=out+"/subset/{sampsize}.psam",
        log=temp(out+"/subset/{sampsize}.log"),
    resources:
        runtime=3,
    log:
        logs + "/subset/{sampsize}",
    benchmark:
        bench + "/subset/{sampsize}",
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --keep <(head -n {params.sampsize} {input.psam} | grep -Ev '^#' | cut -f1) "
        "--make-pgen --pfile {params.prefix} --out {params.out} &>{log}"
