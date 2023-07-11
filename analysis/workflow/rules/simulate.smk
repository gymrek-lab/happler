from pathlib import Path


out = config["out"] + "/simulate"
logs = out + "/logs"
bench = out + "/bench"


locus_chr = config["locus"].split(":")[0]
locus_start = config["locus"].split(":")[1].split('-')[0]
locus_end = config["locus"].split(":")[1].split('-')[1]


rule create_haps_ld_range:
    """ create a hap file with haplotypes that are in a range of LD with their SNPs """
    input:
        gts=Path(config["gts_snp_panel"]).with_suffix(".pvar"),
    output:
        hap=out + "/ld_range/random.hap"
    resources:
        runtime="0:05:00"
    log:
        logs + "/ld_range/create_haps_ld_range",
    benchmark:
        bench + "/ld_range/create_haps_ld_range",
    conda:
        "happler"
    shell:
        "../scripts/choose_different_ld.py {input.gts} > {output.hap}"


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
        hap=out + "/hap/haplotype.hap"
    resources:
        runtime="0:05:00"
    log:
        logs + "/hap/create_hap",
    benchmark:
        bench + "/hap/create_hap",
    conda:
        "../envs/default.yml"
    script:
        "../scripts/create_hap_file.sh"


rule transform:
    """ use the hap file to create a pgen of the haplotype """
    input:
        hap=rules.create_hap.output.hap,
        pgen=config["gts_snp_panel"],
        pvar=Path(config["gts_snp_panel"]).with_suffix(".pvar"),
        psam=Path(config["gts_snp_panel"]).with_suffix(".psam"),
    output:
        pgen=temp(out + "/hap/transform.pgen"),
        pvar=temp(out + "/hap/transform.pvar"),
        psam=temp(out + "/hap/transform.psam"),
    resources:
        runtime="0:05:00"
    log:
        logs + "/hap/transform",
    benchmark:
        bench + "/hap/transform",
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
        pheno=out + "/hap/{beta}.pheno",
    resources:
        runtime="0:05:00"
    log:
        logs + "/hap/{beta}/simphenotype",
    benchmark:
        bench + "/hap/{beta}/simphenotype",
    conda:
        "happler"
    shell:
        "haptools simphenotype -o {output.pheno} {input.pgen} "
        "<( sed 's/\\t0.99$/\\t{params.beta}/' {input.hap} ) &>{log}"
