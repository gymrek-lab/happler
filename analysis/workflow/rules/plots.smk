from pathlib import Path

out = config["out"]
logs = out + "/logs"
bench = out + "/bench"

mode = config["mode"]
out += "/" + mode
logs += "/" + mode
bench += "/" + mode

def agg_ld(wildcards):
    """ a helper function for the other agg functions """
    checkpoint_output = config["ld_range_checkpoint"].get(**wildcards).output.hap
    ld_vals = glob_wildcards(Path(checkpoint_output) / "ld_{ld}/haplotype.hap").ld
    return checkpoint_output, ld_vals

def agg_ld_range_obs(wildcards):
    """ return a list of hap files from happler """
    if "ld_range_checkpoint" in config:
        checkpoint_output, ld_vals = agg_ld(wildcards)
        return expand(
            config["happler_hap"],
            ld=ld_vals,
            alpha=config["mode_attrs"]["alpha"],
            beta=config["mode_attrs"]["beta"],
        )
    else:
        return expand(
            config["happler_hap"],
            beta=config["mode_attrs"]["beta"],
            alpha=config["mode_attrs"]["alpha"],
        )

def agg_ld_range_causal(wildcards):
    """ return a list of hap files from the LD range checkpoint """
    if "ld_range_checkpoint" in config:
        checkpoint_output, ld_vals = agg_ld(wildcards)
        return expand(
            str(checkpoint_output),
            ld=ld_vals,
            beta=config["mode_attrs"]["beta"],
        )
    else:
        return expand(config["causal_hap"], beta=config["mode_attrs"]["beta"])

def agg_ld_range_metrics(wildcards):
    """ return a list of metrics files from the LD range checkpoint """
    if "ld_range_checkpoint" in config:
        checkpoint_output, ld_vals = agg_ld(wildcards)
        return expand(
            config["happler_metrics"],
            ld=ld_vals,
            beta=config["mode_attrs"]["beta"],
            alpha=config["mode_attrs"]["alpha"],
        )
    else:
        return expand(
            config["happler_metrics"],
            beta=config["mode_attrs"]["beta"],
            alpha=config["mode_attrs"]["alpha"],
        )


rule params:
    """ check how wildcards affect the haplotypes output by happler """
    input:
        gts=config["snp_panel"],
        gts_pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        gts_psam=Path(config["snp_panel"]).with_suffix(".psam"),
        observed_haps=agg_ld_range_obs,
        causal_hap=agg_ld_range_causal,
    params:
        observed_haps = lambda wildcards: config["happler_hap"],
        causal_hap = lambda wildcards: config["causal_hap"],
    output:
        png=out + "/happler_params.png",
    resources:
        runtime_min=10,
    log:
        logs + "/plot_params",
    benchmark:
        bench + "/plot_params",
    conda:
        "happler"
    shell:
        "workflow/scripts/parameter_plot.py -o {output.png} "
        "{input.gts} {params.observed_haps} {params.causal_hap} &> {log}"


rule metrics:
    """ check how wildcards affect the finemapped haplotypes output by happler """
    input:
        gts=config["snp_panel"],
        gts_pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        gts_psam=Path(config["snp_panel"]).with_suffix(".psam"),
        observed_haps=agg_ld_range_obs,
        causal_hap=agg_ld_range_causal,
        metrics=agg_ld_range_metrics,
    params:
        observed_haps = lambda wildcards: config["happler_hap"],
        causal_hap = lambda wildcards: config["causal_hap"],
        metrics = lambda wildcards: config["happler_metrics"],
    output:
        png=out + "/finemapping_metrics.png",
    resources:
        runtime_min=10,
    log:
        logs + "/metrics",
    benchmark:
        bench + "/metrics",
    conda:
        "happler"
    shell:
        "workflow/scripts/parameter_plot.py -o {output.png} -m {params.metrics} "
        "{input.gts} {params.observed_haps} {params.causal_hap} &> {log}"
