from pathlib import Path


out = config["out"] + "/simulate"
logs = out + "/logs"
bench = out + "/bench"

mode = config["mode"]
out += "/" + mode
logs += "/" + mode
bench += "/" + mode


locus_chr = config["locus"].split(":")[0]
locus_start = config["locus"].split(":")[1].split('-')[0]
locus_end = config["locus"].split(":")[1].split('-')[1]


checkpoint create_hap_ld_range:
    """ create a hap file suitable for haptools transform and simphenotype """
    input:
        gts=Path(config["gts_snp_panel"]),
        gts_pvar=Path(config["gts_snp_panel"]).with_suffix(".pvar"),
        gts_psam=Path(config["gts_snp_panel"]).with_suffix(".psam"),
    params:
        reps = config["modes"]["ld_range"]["reps"],
        min_ld = config["modes"]["ld_range"]["min_ld"],
        max_ld = config["modes"]["ld_range"]["max_ld"],
        step = config["modes"]["ld_range"]["step"],
        min_af = config["modes"]["ld_range"]["min_af"],
        max_af = config["modes"]["ld_range"]["max_af"],
    output:
        hap=directory(out + "/create_ld_range")
    resources:
        runtime_min=5,
    log:
        logs + "/create_hap_ld_range",
    benchmark:
        bench + "/create_hap_ld_range",
    conda:
        "happler"
    shell:
        "workflow/scripts/choose_different_ld.py -r {params.reps} "
        "--min-ld {params.min_ld} --max-ld {params.max_ld} "
        "--step {params.step} --min-af {params.min_af} "
        "--max-af {params.max_af} {input.gts} {output.hap} &> {log}"


rule create_hap:
    """ create a hap file suitable for haptools transform and simphenotype """
    input:
        gts=Path(config["gts_snp_panel"]).with_suffix(".pvar"),
    params:
        ignore="this", # the first parameter is always ignored for some reason
        chrom=locus_chr,
        locus=config["locus"].split(":")[1].replace('-', '\t'),
        beta=0.99,
        alleles=lambda wildcards: config["modes"]["hap"]["alleles"],
    output:
        hap=out + "/haplotype.hap"
    resources:
        runtime_min=5,
    log:
        logs + "/create_hap",
    benchmark:
        bench + "/create_hap",
    conda:
        "../envs/default.yml"
    script:
        "../scripts/create_hap_file.sh"

transform_input = rules.create_hap.output.hap
if mode == "ld_range":
    transform_input = rules.create_hap_ld_range.output.hap + "/ld_{ld}/haplotype.hap"

if mode == "ld_range":
    out += "/pheno/ld_{ld}"
    logs += "/pheno/ld_{ld}"
    bench += "/pheno/ld_{ld}"

rule transform:
    """ use the hap file to create a pgen of the haplotype """
    input:
        hap=transform_input,
        pgen=config["gts_snp_panel"],
        pvar=Path(config["gts_snp_panel"]).with_suffix(".pvar"),
        psam=Path(config["gts_snp_panel"]).with_suffix(".psam"),
    output:
        pgen=temp(out + "/transform.pgen"),
        pvar=temp(out + "/transform.pvar"),
        psam=temp(out + "/transform.psam"),
    resources:
        runtime_min=5,
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
    output:
        pheno=out + "/{beta}.pheno",
    resources:
        runtime_min=5,
    log:
        logs + "/{beta}/simphenotype",
    benchmark:
        bench + "/{beta}/simphenotype",
    conda:
        "happler"
    shell:
        "haptools simphenotype -o {output.pheno} {input.pgen} "
        "<( sed 's/\\t0.99$/\\t{params.beta}/' {input.hap} ) &>{log}"
