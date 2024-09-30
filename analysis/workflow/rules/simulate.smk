from pathlib import Path


out = config["out"] + "/simulate"
logs = out + "/logs"
bench = out + "/bench"

mode = config["mode"]
out += "/" + mode
logs += "/" + mode
bench += "/" + mode

if mode == "ld_range" and config["random"]:
    out += "/random"
    logs += "/random"
    bench += "/random"

def parse_locus(locus):
    """parse locus into chrom, start, end"""
    chrom = locus.split("_")[0]
    end = locus.split("-")[1]
    start = locus.split("_")[1].split("-")[0]
    return chrom, start, end

hap_ld_range_output = out + "/create_ld_range/ld_{ld}/haplotype.hap"

checkpoint create_hap_ld_range:
    """ create a hap file suitable for haptools transform and simphenotype """
    input:
        gts=Path(config["gts_snp_panel"]),
        gts_pvar=Path(config["gts_snp_panel"]).with_suffix(".pvar"),
        gts_psam=Path(config["gts_snp_panel"]).with_suffix(".psam"),
    params:
        min_ld = config["modes"]["ld_range"]["min_ld"],
        max_ld = config["modes"]["ld_range"]["max_ld"],
        step = config["modes"]["ld_range"]["step"],
        min_af = config["modes"]["ld_range"]["min_af"],
        max_af = config["modes"]["ld_range"]["max_af"],
        out = lambda wildcards: expand(
            hap_ld_range_output, **wildcards, allow_missing=True,
        ),
        seed = 42,
    output:
        hap=directory(out + "/create_ld_range")
    resources:
        runtime=5,
    log:
        logs + "/create_hap_ld_range",
    benchmark:
        bench + "/create_hap_ld_range",
    conda:
        "happler"
    shell:
        "workflow/scripts/choose_different_ld.py "
        "--min-ld {params.min_ld} --max-ld {params.max_ld} "
        "--step {params.step} --min-af {params.min_af} --seed {params.seed} "
        "--max-af {params.max_af} {input.gts} {params.out} &> {log}"


rule create_hap:
    """ create a hap file suitable for haptools transform and simphenotype """
    input:
        gts=Path(config["gts_snp_panel"]).with_suffix(".pvar"),
    params:
        ignore="this", # the first parameter is always ignored for some reason
        chrom=lambda wildcards: parse_locus(wildcards.locus)[0],
        locus=lambda wildcards: wildcards.locus.split("_")[1].replace('-', '\t'),
        beta=0.99,
        alleles=lambda wildcards: config["modes"]["hap"]["alleles"],
    output:
        hap=out + "/haplotype.hap"
    resources:
        runtime=5,
    log:
        logs + "/create_hap",
    benchmark:
        bench + "/create_hap",
    conda:
        "../envs/default.yml"
    script:
        "../scripts/create_hap_file.sh"

if mode == "ld_range":
    out += "/pheno/ld_{ld}"
    logs += "/pheno/ld_{ld}"
    bench += "/pheno/ld_{ld}"

rule transform:
    """ use the hap file to create a pgen of the haplotype """
    input:
        hap=hap_ld_range_output if mode == "ld_range" else rules.create_hap.output.hap,
        pgen=config["gts_snp_panel"],
        pvar=Path(config["gts_snp_panel"]).with_suffix(".pvar"),
        psam=Path(config["gts_snp_panel"]).with_suffix(".psam"),
    output:
        pgen=out + "/transform.pgen",
        pvar=out + "/transform.pvar",
        psam=out + "/transform.psam",
    resources:
        runtime=5,
    log:
        logs + "/transform",
    benchmark:
        bench + "/transform",
    conda:
        "happler"
    shell:
        "haptools transform -o {output.pgen} {input.pgen} {input.hap} &>{log}"


rule simphenotype:
    """ use the hap file to create simulated phenotypes for the haplotype """
    input:
        hap=rules.transform.input.hap,
        pgen=rules.transform.output.pgen,
        pvar=rules.transform.output.pvar,
        psam=rules.transform.output.psam,
    params:
        beta=lambda wildcards: wildcards.beta,
        seed = 42,
    output:
        pheno=out + "/{beta}.pheno",
    resources:
        runtime=5,
    log:
        logs + "/{beta}/simphenotype",
    benchmark:
        bench + "/{beta}/simphenotype",
    conda:
        "happler"
    shell:
        "haptools simphenotype --seed {params.seed} -o {output.pheno} {input.pgen} "
        "<( sed 's/\\t0.99$/\\t{params.beta}/' {input.hap} ) &>{log}"
