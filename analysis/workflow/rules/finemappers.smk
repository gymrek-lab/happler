from pathlib import Path


out = config["out"]
logs = out + "/logs"
bench = out + "/bench"


exclude_causal = {"in": 0, "ex": 1}


rule susie:
    """run the method SuSiE"""
    input:
        gt=config["snp_panel"],
        phen=config["pheno"],
    params:
        outdir=lambda wildcards, output: Path(output.susie).parent,
        exclude_causal=lambda wildcards: expand(config["causal_gt"], **wildcards)[0] if \
            int(exclude_causal[wildcards.causal]) else "NULL",
    output:
        susie=out + "/{causal}clude/susie.rds",
    resources:
        runtime=15,
        queue="hotel",
    log:
        logs + "/{causal}clude/susie",
    benchmark:
        bench + "/{causal}clude/susie",
    conda:
        "../envs/susie.yml"
    shell:
        "workflow/scripts/run_SuSiE.R {input} {params} &>{log}"


rule finemap:
    """run the method FINEMAP"""
    input:
        gt=config["snp_panel"],
        phen=config["pheno"],
    params:
        outdir=lambda wildcards, output: Path(output.finemap).parent,
        exclude_causal=lambda wildcards: expand(config["causal_gt"], **wildcards)[0] if \
            not int(exclude_causal[wildcards.causal]) else "NULL",
    output:
        sumstats=temp(out + "/{causal}clude/sumstats.rds"),
        finemap=out + "/{causal}clude/finemap.rds",
    resources:
        runtime=45,
        queue="hotel",
    log:
        logs + "/{causal}clude/finemap",
    benchmark:
        bench + "/{causal}clude/finemap",
    conda:
        "../envs/susie.yml"
    shell:
        "workflow/scripts/run_FINEMAP.R {input} {params} &>{log}"


rule results:
    """create plots to summarize the results of the simulations"""
    input:
        gt=config["snp_panel"],
        finemap=rules.finemap.output.finemap,
        susie=rules.susie.output.susie,
        causal_gt=config["causal_gt"],
    params:
        outdir=lambda wildcards, output: Path(output.susie_pdf).parent,
        causal_hap="",
    output:
        finemap_pdf= out + "/{causal}clude/finemap.pdf",
        susie_pdf= out + "/{causal}clude/susie.pdf",
    resources:
        runtime=60,
        queue="hotel",
        # susie_eff_pdf=temp( out + "/{causal}clude/susie_eff.pdf"),
    log:
        logs + "/{causal}clude/results",
    conda:
        "../envs/susie.yml"
    script:
        "../scripts/summarize_results.R"
