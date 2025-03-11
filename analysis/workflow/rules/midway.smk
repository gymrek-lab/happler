from pathlib import Path

out = config["out"]
logs = out + "/logs"
bench = out + "/bench"


rule sub_pheno:
    """
    subset the phenotype file to include only the desired replicate
    This rule is copied from happler.smk
    """
    input:
        pheno=config["pheno"],
    params:
        rep=lambda wildcards: int(wildcards.rep)+2,
    output:
        pheno=out + "/phen.pheno",
    wildcard_constraints:
        rep="\d+"
    resources:
        runtime=7,
    threads: 1,
    log:
        logs + "/sub_pheno",
    benchmark:
        bench + "/sub_pheno",
    conda:
        "../envs/default.yml"
    shell:
        "cut -f 1,{params.rep} {input.pheno} | "
        "(echo -e \"#IID\\thap\" && tail -n+2) >{output.pheno} 2>{log}"


pheno = rules.sub_pheno.output.pheno if "{rep}" in out else config["pheno"]


rule manhattan:
    """ run happler midway-through via a plink2 GWAS """
    input:
        gts=config["snp_panel"],
        gts_pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        gts_psam=Path(config["snp_panel"]).with_suffix(".psam"),
        pts=pheno,
        hap=config["hap_file"]
    params:
        out_prefix = lambda wildcards, output: str(output.dir) + "/out",
        maf=config["min_maf"],
        target = "H0",
        tswitch=lambda wildcards: 2 if wildcards.switch == "tscore" else 1,
    output:
        dir=directory(out + "/{switch}"),
        linear=out + "/{switch}/out.linear",
    wildcard_constraints:
        switch="(tscore|interact)"
    resources:
        runtime=10,
    log:
        logs + "/{switch}/manhattan",
    benchmark:
        bench + "/{switch}/manhattan",
    conda:
        "happler"
    shell:
        "rsid=\"$(grep -E '^V' {input.hap} | cut -f5 | tail -n1)\" && "
        "workflow/scripts/midway_manhattan.bash {input.gts} {input.pts} {input.hap} "
        "{params.out_prefix} {params.target} \"$rsid\" {params.maf} {params.tswitch} "
        "&> {log}"
