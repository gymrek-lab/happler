from pathlib import Path


out = config["out"] + "/finemappers"
logs = out + "/logs"
bench = out + "/bench"


exclude_causal = {"in": 0, "ex": 1}


rule run:
    """run the methods FINEMAP and SuSie"""
    input:
        gt=config["snp_panel"],
        phen=config["pheno"],
    params:
        outdir=lambda wildcards, output: Path(output.susie).parent,
        exclude_causal=lambda wildcards: int(exclude_causal[wildcards.causal]),
    output:
        sumstats=temp(out + "/{causal}clude/sumstats.rds"),
        finemap=out + "/{causal}clude/finemap.rds",
        susie=out + "/{causal}clude/susie.rds",
    resources:
        runtime="1:00:00",
    log:
        logs + "/{causal}clude/run",
    benchmark:
        bench + "/{causal}clude/run",
    conda:
        "../envs/susie.yml"
    shell:
        "workflow/scripts/finemapping_methods.R {input} {params} &>{log}"


rule results:
    """create plots to summarize the results of the simulations"""
    input:
        gt=rules.run.input.gt,
        finemap=rules.run.output.finemap,
        susie=rules.run.output.susie,
    params:
        outdir=lambda wildcards, output: Path(output.susie_pdf).parent,
        exclude_causal=lambda wildcards: int(exclude_causal[wildcards.causal]),
        causal_hap="",
    output:
        finemap_pdf= out + "/{causal}clude/finemap.pdf",
        susie_pdf= out + "/{causal}clude/susie.pdf",
        # susie_eff_pdf=temp( out + "/{causal}clude/susie_eff.pdf"),
    log:
        logs + "/{causal}clude/results",
    conda:
        "../envs/susie.yml"
    script:
        "scripts/summarize_results.R"
