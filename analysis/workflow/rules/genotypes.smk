from pathlib import Path


out = "out/genotypes"
logs = "out/logs/genotypes"
bench = "out/bench/genotypes"


locus_chr = config["locus"].split(":")[0]
locus_start = config["locus"].split("-")[1].split('-')[0]
locus_end = config["locus"].split("-")[1].split('-')[1]



rule plink2vcf:
    """ convert a PLINK file to VCF """
    input:
        pgen=config["snp_panel"],
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
        runtime="0:15:00"
    log:
        logs + "/plink2vcf"
    benchmark:
        bench + "/plink2vcf"
    conda:
        "envs/default.yml"
    shell:
        "plink2 --pfile {params.pfile} --out {params.out} --from-bp {params.start} "
        "--to-bp {params.end} --chr {params.chrom} --snps-only 'just-acgt' "
        "--export vcf bgz id-paste=iid &>{log} && "
        "bcftools view -Ob -o {output.bcf} {output.vcf} &>>{log} && "
        "tabix -p bcf {output.bcf} &>>{log}"


def phase_gt_input(wildcards):
    input_files = { 'map': config['phase_map'] }
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
        "envs/shapeit.yml"
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
        str_vcf=config["str_panel"],
        str_vcf_idx=config["str_panel"],
        samp=config["exclude_samples"],
    output:
        samples=out+"/samples.tsv"
    log:
        logs+"/keep_samps"
    conda:
        "envs/default.yml"
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
        samples=rules.keep_samps.output.samples,
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
        "envs/default.yml"
    shell:
        "plink2 --bcf {input.vcf} --maf {params.maf} --make-pgen "
        "--keep {input.samples} --out {params.prefix} &>{log}"
