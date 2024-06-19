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
        out=lambda wildcards, output: str(Path(output.log).with_suffix('')),
        start=lambda wildcards: parse_locus(wildcards.locus)[1],
        end=lambda wildcards: parse_locus(wildcards.locus)[2],
        chrom=lambda wildcards: parse_locus(wildcards.locus)[0],
    output:
        vcf=temp(out + "/unphased.vcf.gz"),
        vcf_idx=temp(out + "/unphased.vcf.gz.tbi"),
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
        "tabix -p vcf {output.vcf} &>>{log}"


def phase_gt_input(wildcards):
    chrom = parse_locus(wildcards.locus)[0]
    input_files = { 'map': config['phase_map'].format(chr=chrom) }
    if config["snp_panel"].endswith('.pgen'):
        input_files['unphased'] = rules.plink2vcf.output.vcf
    else:
        input_files['unphased'] = config["snp_panel"].format(chr=chrom)
    input_files['unphased_idx'] = input_files['unphased'] + ".tbi"
    input_files["ref"] = config["ref_panel"].format(chr=chrom)
    input_files["log"] = rules.plink2vcf.output.log
    return input_files


def get_num_variants(wildcards):
    with open(expand(rules.plink2vcf.output.log, **wildcards)[0], "r") as file:
        return int(next(line for line in file if 'variants rem' in line).split(" ")[0])


rule phase_gt:
    """ phase an unphased set of genotypes """
    input: unpack(phase_gt_input)
    params:
        locus=lambda wildcards: wildcards.locus.replace("_", ":"),
        prefix=lambda wildcards, output: str(Path(output.log).with_suffix("")),
    output:
        phased=out + "/phased.vcf.gz",
        phased_idx=out + "/phased.vcf.gz.tbi",
        log=temp(out + "/phased.log"),
    resources:
        # We use a custom formula to determine the time requirements:
        # This computes the number of minutes from the number of variants
        runtime=lambda wildcards: int(
            0.0015600416766108 * get_num_variants(wildcards) + 57.27754932492151
        ),
        # We use a custom formula to determine the memory requirements:
        # This computes the number of GB from the number of variants and then it
        # multiplies by 1000 MB per GB
        mem_mb=lambda wildcards: int(
            (0.0006264258500361132 * get_num_variants(wildcards) + 13.0220968394949) * 1000
        ),
    threads: 32
    log:
        logs + "/phase_gt"
    benchmark:
        bench + "/phase_gt"
    conda:
        "../envs/beagle.yml"
    shell:
        "beagle -Xmx{resources.mem_mb}m gt={input.unphased} ref={input.ref} "
        "out={params.prefix} map={input.map} chrom={params.locus} nthreads={threads} "
        "impute=false &>{log} && "
        "tabix -p vcf {output.phased} &>>{log}"


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
        samples=(
            (rules.keep_samps.output.samples if check_config("str_panel") else config["exclude_samples"])
            if check_config('exclude_samples') else []
        ),
    params:
        maf=config["min_maf"],
        prefix=lambda wildcards, output: Path(output.pgen).with_suffix(""),
        samps=lambda wildcards, input: (" --" + (
            "keep " if check_config("str_panel") else "remove "
        ) + input.samples) if check_config("exclude_samples") else "",
    output:
        pgen=out+"/snps.pgen",
        pvar=out+"/snps.pvar",
        psam=out+"/snps.psam",
        log=temp(out+"/snps.log"),
    resources:
        runtime=15,
    threads: 2
    log:
        logs + "/vcf2plink",
    benchmark:
        bench + "/vcf2plink",
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --vcf {input.vcf} --maf {params.maf} --geno 0 --make-pgen "
        "--threads {threads}{params.samps} --out {params.prefix} &>{log}"


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
    output:
        pgen=out+"/subset/{sampsize}.pgen",
        pvar=out+"/subset/{sampsize}.pvar",
        psam=out+"/subset/{sampsize}.psam",
        log=temp(out+"/subset/{sampsize}.log"),
    resources:
        runtime=12,
    threads: 2
    log:
        logs + "/subset/{sampsize}",
    benchmark:
        bench + "/subset/{sampsize}",
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --keep <(head -n {params.sampsize} {input.psam} | grep -Ev '^#' | cut -f1) "
        "--make-pgen --pfile {params.prefix} --out {params.out} &>{log}"
