import re
import sys
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
    ld_vals = glob_wildcards(Path(checkpoint_output) / "{num_haps}_haps/ld_{ld}/haplotype.hap").ld
    return checkpoint_output, sorted(list(set(ld_vals)))

def agg_ld_range_obs(wildcards):
    """ return a list of hap files from happler """
    if "ld_range_checkpoint" in config:
        checkpoint_output, ld_vals = agg_ld(wildcards)
        return expand(
            config["happler_hap"],
            ld=ld_vals,
            alpha=config["mode_attrs"]["alpha"],
            beta=config["mode_attrs"]["beta"],
            num_haps=config["mode_attrs"]["num_haps"],
            rep=range(config["mode_attrs"]["reps"]),
            **wildcards,
        )
    else:
        return expand(
            config["happler_hap"],
            beta=config["mode_attrs"]["beta"],
            num_haps=config["mode_attrs"]["num_haps"],
            **wildcards,
        )

def agg_ld_range_causal(wildcards):
    """ return a list of hap files from the LD range checkpoint """
    if "ld_range_checkpoint" in config:
        checkpoint_output, ld_vals = agg_ld(wildcards)
        return expand(
            str(checkpoint_output),
            ld=ld_vals,
            beta=config["mode_attrs"]["beta"],
            num_haps=config["mode_attrs"]["num_haps"],
            **wildcards,
        )
    elif mode == "midway":
        return expand(
            config["causal_hap"],
            locus = config["loci"],
        )
    else:
        return expand(
            config["causal_hap"],
            beta=config["mode_attrs"]["beta"],
            num_haps=config["mode_attrs"]["num_haps"],
            **wildcards,
        )

def agg_ld_range_metrics(wildcards):
    """ return a list of metrics files from the LD range checkpoint """
    if "ld_range_checkpoint" in config:
        checkpoint_output, ld_vals = agg_ld(wildcards)
        return expand(
            config["happler_metrics"],
            ld=ld_vals,
            beta=config["mode_attrs"]["beta"],
            num_haps=config["mode_attrs"]["num_haps"],
            alpha=config["mode_attrs"]["alpha"],
            rep=range(config["mode_attrs"]["reps"]),
            **wildcards,
        )
    else:
        return expand(
            config["happler_metrics"],
            beta=config["mode_attrs"]["beta"],
            num_haps=config["mode_attrs"]["num_haps"],
            **wildcards,
        )


switch_sim_mode = {
    "interact": ("hap", "indep"),
    "tscore": ("hap", "parent"),
}


def agg_midway_linear(wildcards):
    """ return a list of midway linear files """
    sim_modes = switch_sim_mode[wildcards.switch]
    return expand(
        config["midway_linear"],
        locus=config["loci"],
        rep=range(config["mode_attrs"]["reps"]),
        sim_mode=sim_modes,
        beta=config["mode_attrs"]["beta"],
        **wildcards,
        allow_missing=True,
    )

fill_out_globals = lambda wildcards, val: expand(
    val,
    locus=wildcards.locus,
    sampsize=wildcards.sampsize,
    allow_missing=True,
)

fill_out_globals_midway = lambda wildcards, val: expand(
    val,
    sampsize=wildcards.sampsize,
    switch=wildcards.switch,
    allow_missing=True,
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
        observed_haps = lambda wildcards: fill_out_globals(wildcards, config["happler_hap"]),
        causal_hap = lambda wildcards: fill_out_globals(wildcards, config["causal_hap"]),
    output:
        png=out + "/happler_params.png",
    resources:
        runtime=20,
    log:
        logs + "/plot_params",
    benchmark:
        bench + "/plot_params",
    conda:
        "happler"
    shell:
        "workflow/scripts/parameter_plot.py -o {output.png} "
        "--order num_haps,beta,ld "
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
        observed_haps = lambda wildcards: fill_out_globals(wildcards, config["happler_hap"]),
        causal_hap = lambda wildcards: fill_out_globals(wildcards, config["causal_hap"]),
        metrics = lambda wildcards: fill_out_globals(wildcards, config["happler_metrics"]),
    output:
        png=out + "/finemapping_metrics.png",
    resources:
        runtime=20,
    log:
        logs + "/metrics",
    benchmark:
        bench + "/metrics",
    conda:
        "happler"
    shell:
        "workflow/scripts/parameter_plot.py -o {output.png} -m {params.metrics} "
        "--order num_haps,beta,ld "
        "{input.gts} {params.observed_haps} {params.causal_hap} &> {log}"


create_glob_from_wildcards = lambda path, wildcards: re.sub(r"\{[^}]+\}", "*", expand(path, sim_mode="{"+switch_sim_mode[wildcards.switch].join(",")+"}", allow_missing=True))


rule midway:
    """summarize the results from many midway-manhattan runs"""
    input:
        linears=agg_midway_linear,
        snplists=agg_ld_range_causal,
    params:
        case_type="sim_mode",
        pos_type="hap",
        linears=lambda wildcards: fill_out_globals_midway(wildcards, config["midway_linear"]),
        causal_hap = lambda wildcards: fill_out_globals_midway(wildcards, config["causal_hap"]),
        linears_glob = lambda wildcards: expand(
            re.sub(
                r"\{(?!sim_mode\})[^}]+\}", "*",
                fill_out_globals_midway(wildcards, config["midway_linear"])[0],
            ),
            sim_mode="{"+",".join(switch_sim_mode[wildcards.switch])+"}",
            allow_missing=True,
        )[0]
    output:
        png=out + "/midway_summary.{switch}.png",
    resources:
        runtime=7,
    log:
        logs + "/midway.{switch}",
    benchmark:
        bench + "/midway.{switch}",
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/midway_manhattan_summary.py "
        "-o {output.png} --verbosity DEBUG --pos-type {params.pos_type} "
        "-f <(ls -1 {params.linears_glob}) --color locus "
        "{params.linears} {params.causal_hap} {params.case_type} &>{log}"
