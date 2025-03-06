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
elif mode == "midway":
    out += "/sim_{sim_mode}"
    logs += "/sim_{sim_mode}"
    bench += "/sim_{sim_mode}"

def parse_locus(locus):
    """parse locus into chrom, start, end"""
    chrom = locus.split("_")[0]
    end = locus.split("-")[1]
    start = locus.split("_")[1].split("-")[0]
    return chrom, start, end

hap_ld_range_output = out + "/create_ld_range/{num_haps}_haps/ld_{ld}/haplotype.hap"

gts_panel = Path(
    config["gts_str_panel"] if mode == "str" else config["gts_snp_panel"]
)

checkpoint create_hap_ld_range:
    """ create a hap file suitable for haptools transform and simphenotype """
    input:
        gts=Path(config["gts_snp_panel"]),
        gts_pvar=Path(config["gts_snp_panel"]).with_suffix(".pvar"),
        gts_psam=Path(config["gts_snp_panel"]).with_suffix(".psam"),
    params:
        min_ld = lambda wildcards: config["modes"]["ld_range"]["min_ld"],
        max_ld = lambda wildcards: config["modes"]["ld_range"]["max_ld"],
        step = lambda wildcards: config["modes"]["ld_range"]["step"],
        min_af = lambda wildcards: config["modes"]["ld_range"]["min_af"],
        max_af = lambda wildcards: config["modes"]["ld_range"]["max_af"],
        num_haps = lambda wildcards: "-n " + " -n ".join(map(str, config["modes"]["ld_range"]["num_haps"])),
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
        "workflow/scripts/choose_different_ld.py {params.num_haps} "
        "--min-ld {params.min_ld} --max-ld {params.max_ld} "
        "--step {params.step} --min-af {params.min_af} --seed {params.seed} "
        "--max-af {params.max_af} {input.gts} {params.out} &> {log}"


rule create_hap:
    """ create a hap file suitable for haptools transform and simphenotype """
    input:
        gts=gts_panel.with_suffix(".pvar"),
    params:
        ignore="this", # the first parameter is always ignored for some reason
        chrom=lambda wildcards: parse_locus(wildcards.locus)[0],
        locus=lambda wildcards: wildcards.locus.split("_")[1].replace('-', '\t'),
        beta=0.99,
        alleles=lambda wildcards: config["modes"]["hap"]["alleles"] if mode == "hap" else [],
        repeat=lambda wildcards: config["modes"]["str"]["id"] if mode == "str" else 0,
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

transform_input = rules.create_hap.output.hap
if mode == "ld_range":
    out += "/pheno/{num_haps}_haps/ld_{ld}"
    logs += "/pheno/{num_haps}_haps/ld_{ld}"
    bench += "/pheno/{num_haps}_haps/ld_{ld}"
    transform_input = hap_ld_range_output
elif mode == "midway":
    transform_input = lambda wildcards: expand(
        config["modes"][mode]["haps"] + "/{fname}", fname=config["locus_traits"][wildcards.locus],
    )[0]

rule transform:
    """ use the hap file to create a pgen of the haplotype """
    input:
        hap=transform_input,
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


make_midway_hap = lambda sim_mode: "sed '${/^V/ d;}'" if sim_mode == "wo" else "cat"


rule simphenotype:
    """ use the hap file to create simulated phenotypes for the haplotype """
    input:
        hap=rules.transform.input.hap,
        pgen=rules.transform.output.pgen if mode != "str" else config["gts_str_panel"],
        pvar=rules.transform.output.pvar if mode != "str" else Path(config["gts_str_panel"]).with_suffix(".pvar"),
        psam=rules.transform.output.psam if mode != "str" else Path(config["gts_str_panel"]).with_suffix(".psam"),
    params:
        # if the beta value is 0, we just use the beta of the hap file
        # otherwise, we change it on the fly
        hap=lambda wildcards, input: "<(" + make_midway_hap(wildcards.sim_mode) + " " + input.hap + (
            "" if wildcards.beta == "0" else
            " | awk -F $'\\t' -v 'OFS=\\t' \'$1 == \"H\" { $6 = " + wildcards.beta + " }1\'"
        )+")",
        seed = 42,
        reps = lambda wildcards: config["modes"]["ld_range"]["reps"],
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
        "haptools simphenotype --seed {params.seed} -o {output.pheno} "
        "--replications {params.reps} {input.pgen} {params.hap} &>{log}"
