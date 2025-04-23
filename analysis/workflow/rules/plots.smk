import re
import sys
from pathlib import Path
from functools import partial

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
    "covariance": ("hap", "parent"),
    "bic": ("hap", "parent"),
    "interact-bic": ("hap", "indep"),
    "pip-parent": ("hap", "parent"),
    "pip-interact": ("hap", "indep"),
}


def agg_midway_linear(wildcards, beta: bool = False):
    """ return a list of midway linear files """
    sim_modes = switch_sim_mode[wildcards.switch]
    expand_partial = expand
    if beta:
        expand_partial = partial(expand_partial, beta=config["mode_attrs"]["beta"])
    else:
        expand_partial = partial(expand_partial, beta=wildcards.beta)
    return expand_partial(
        config["midway_linear"],
        locus=config["loci"],
        rep=range(config["mode_attrs"]["reps"]),
        sim_mode=sim_modes,
        **wildcards,
        allow_missing=True,
    )


def agg_finemap(wildcards, beta: bool = False, also_exclude = False):
    """ return a list of finemap files """
    sim_modes = switch_sim_mode[wildcards.switch]
    expand_partial = expand
    if beta:
        expand_partial = partial(expand_partial, beta=config["mode_attrs"]["beta"])
    else:
        expand_partial = partial(expand_partial, beta=wildcards.beta)
    if also_exclude:
        return expand_partial(
            config["finemap_metrics"],
            locus=config["loci"],
            rep=range(config["mode_attrs"]["reps"]),
            sim_mode=("hap",),
            ex=("ex","in",),
            **wildcards,
            allow_missing=True,
        )
    else:
        return expand_partial(
            config["finemap_metrics"],
            locus=config["loci"],
            rep=range(config["mode_attrs"]["reps"]),
            sim_mode=sim_modes,
            ex=("in",),
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

fill_out_globals_finemap = lambda wildcards, val: expand(
    val,
    sampsize=wildcards.sampsize,
    switch=wildcards.switch,
    ex=("in",),
    allow_missing=True,
)

fill_out_globals_finemap_beta = lambda wildcards, val: expand(
    fill_out_globals_finemap(wildcards, val),
    beta=wildcards.beta,
    allow_missing=True,
)

fill_out_globals_finemap_cs_len = lambda wildcards, val: expand(
    val,
    sampsize=wildcards.sampsize,
    switch=wildcards.switch,
    sim_mode=("hap",),
    allow_missing=True,
)

fill_out_globals_midway_beta = lambda wildcards, val: expand(
    fill_out_globals_midway(wildcards, val),
    beta=wildcards.beta,
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
        linears=partial(agg_midway_linear, beta=True),
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
        )[0],
        bic=lambda wildcards: "--kind bic " if wildcards.switch in ("bic", "interact-bic") else "",
        thresh=lambda wildcards: "--thresh 3 " if wildcards.switch in ("bic", "interact-bic") else "--thresh 0.05 ",
    output:
        png=out + "/midway_summary.{switch}.pdf",
        metrics=out+"/midway_summary_metrics.{switch}.tsv",
    resources:
        runtime=7,
    wildcard_constraints:
        switch="(interact|tscore|covariance|bic|interact-bic)"
    log:
        logs + "/midway.{switch}",
    benchmark:
        bench + "/midway.{switch}",
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/midway_manhattan_summary.py {params.bic}"
        "-o {output.png} --verbosity DEBUG --pos-type {params.pos_type} "
        "-f <(ls -1 {params.linears_glob}) --color locus {params.thresh}"
        "{params.linears} {params.causal_hap} {params.case_type} >{output.metrics} 2>{log}"


rule midway_beta:
    """summarize the results from many midway-manhattan runs for a specific value of beta"""
    input:
        linears=agg_midway_linear,
        snplists=agg_ld_range_causal,
    params:
        case_type="sim_mode",
        pos_type="hap",
        linears=lambda wildcards: fill_out_globals_midway_beta(wildcards, config["midway_linear"]),
        causal_hap = lambda wildcards: fill_out_globals_midway_beta(wildcards, config["causal_hap"]),
        linears_glob = lambda wildcards: expand(
            re.sub(
                r"\{(?!sim_mode\})[^}]+\}", "*",
                fill_out_globals_midway_beta(wildcards, config["midway_linear"])[0],
            ),
            sim_mode="{"+",".join(switch_sim_mode[wildcards.switch])+"}",
            allow_missing=True,
        )[0],
        bic=lambda wildcards: "--kind bic " if wildcards.switch in ("bic", "interact-bic") else "",
        thresh=lambda wildcards: "--thresh 3 " if wildcards.switch in ("bic", "interact-bic") else "--thresh 0.05 ",
    output:
        png=out + "/beta_{beta}/midway_summary.{switch}.pdf",
        metrics=out+"/beta_{beta}/midway_summary_metrics.{switch}.tsv",
    resources:
        runtime=7,
    wildcard_constraints:
        switch="(interact|tscore|covariance|bic|interact-bic)"
    log:
        logs + "/beta_{beta}/midway.{switch}",
    benchmark:
        bench + "/beta_{beta}/midway.{switch}",
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/midway_manhattan_summary.py {params.bic}"
        "-o {output.png} --verbosity DEBUG --pos-type {params.pos_type} "
        "-f <(ls -1 {params.linears_glob}) --color locus {params.thresh}"
        "{params.linears} {params.causal_hap} {params.case_type} >{output.metrics} 2>{log}"


rule finemap:
    """summarize the results from many finemapping runs"""
    input:
        finemaps=partial(agg_finemap, beta=True),
        haps=agg_ld_range_causal,
    params:
        case_type="sim_mode",
        pos_type="hap",
        finemaps=lambda wildcards: fill_out_globals_finemap(wildcards, config["finemap_metrics"]),
        causal_hap = lambda wildcards: fill_out_globals_midway(wildcards, config["causal_hap"]),
        finemaps_glob = lambda wildcards: expand(
            re.sub(
                r"\{(?!sim_mode\})[^}]+\}", "*",
                fill_out_globals_finemap(wildcards, config["finemap_metrics"])[0],
            ),
            sim_mode="{"+",".join(switch_sim_mode[wildcards.switch])+"}",
            allow_missing=True,
        )[0],
        thresh=lambda wildcards: "--thresh 0.9 ",
    output:
        png=out + "/midway_summary.{switch}.pdf",
        metrics=out+"/midway_summary_metrics.{switch}.tsv",
    resources:
        runtime=7,
    wildcard_constraints:
        switch="(pip-parent|pip-interact)"
    log:
        logs + "/midway.{switch}",
    benchmark:
        bench + "/midway.{switch}",
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/midway_manhattan_summary.py --kind pip {params.thresh}"
        "-o {output.png} --verbosity DEBUG --pos-type {params.pos_type} "
        "-f <(ls -1 {params.finemaps_glob}) --color locus "
        "{params.finemaps} {params.causal_hap} {params.case_type} >{output.metrics} 2>{log}"


rule finemap_beta:
    """summarize the results from many finemapping runs for a specific value of beta"""
    input:
        finemaps=partial(agg_finemap),
        haps=agg_ld_range_causal,
    params:
        case_type="sim_mode",
        pos_type="hap",
        finemaps=lambda wildcards: fill_out_globals_finemap_beta(wildcards, config["finemap_metrics"]),
        causal_hap = lambda wildcards: fill_out_globals_midway_beta(wildcards, config["causal_hap"]),
        finemaps_glob = lambda wildcards: expand(
            re.sub(
                r"\{(?!sim_mode\})[^}]+\}", "*",
                fill_out_globals_finemap_beta(wildcards, config["finemap_metrics"])[0],
            ),
            sim_mode="{"+",".join(switch_sim_mode[wildcards.switch])+"}",
            allow_missing=True,
        )[0],
        thresh=lambda wildcards: "--thresh 0.9 ",
    output:
        png=out + "/beta_{beta}/midway_summary.{switch}.pdf",
        metrics=out+"/beta_{beta}/midway_summary_metrics.{switch}.tsv",
    resources:
        runtime=7,
    wildcard_constraints:
        switch="(pip-parent|pip-interact)"
    log:
        logs + "/beta_{beta}/midway.{switch}",
    benchmark:
        bench + "/beta_{beta}/midway.{switch}",
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/midway_manhattan_summary.py --kind pip {params.thresh}"
        "-o {output.png} --verbosity DEBUG --pos-type {params.pos_type} "
        "-f <(ls -1 {params.finemaps_glob}) --color locus "
        "{params.finemaps} {params.causal_hap} {params.case_type} >{output.metrics} 2>{log}"


rule finemap_cs_length:
    """summarize the results from many finemapping runs. Show CS length instead of PIPs"""
    input:
        finemaps=partial(agg_finemap, beta=True, also_exclude=True),
        haps=agg_ld_range_causal,
    params:
        case_type="ex",
        finemaps=lambda wildcards: fill_out_globals_finemap_cs_len(wildcards, config["finemap_metrics"]),
        causal_hap = lambda wildcards: fill_out_globals_midway(wildcards, config["causal_hap"]),
        finemaps_glob = lambda wildcards: expand(
            re.sub(
                r"\{(?!ex\})[^}]+\}", "*",
                fill_out_globals_finemap_cs_len(wildcards, config["finemap_metrics"])[0],
            ),
            ex="{"+",".join(("ex", "in"))+"}",
            allow_missing=True,
        )[0],
    output:
        png=out + "/cs_len_summary.{switch}.pdf",
        metrics=out+"/cs_len_summary_metrics.{switch}.tsv",
    resources:
        runtime=7,
    wildcard_constraints:
        switch="pip-parent"
    log:
        logs + "/cs_len.{switch}",
    benchmark:
        bench + "/cs_len.{switch}",
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/midway_manhattan_summary.py --kind cs_length "
        "-o {output.png} --verbosity DEBUG "
        "-f <(ls -1 {params.finemaps_glob}) --color locus "
        "{params.finemaps} {params.causal_hap} {params.case_type} >{output.metrics} 2>{log}"


out = out.replace("{sampsize}samples/", "")
logs = logs.replace("{sampsize}samples/", "")
bench = bench.replace("{sampsize}samples/", "")


rule midway_metrics:
    """summarize metrics from many executions of midway_beta"""
    input:
        metrics = expand(
            rules.midway_beta.output.metrics,
            beta=config["mode_attrs"]["beta"],
            sampsize=config["sample_size"],
            allow_missing=True
        ),
    params:
        metrics = lambda wildcards: expand(rules.midway_beta.output.metrics, switch=wildcards.switch, allow_missing=True),
        use_flex_axes = lambda wildcards: "--use-flex-axes-limits " if wildcards.switch == "interact" else "",
    output:
        png=out + "/midway_summary_metrics.{switch}.pdf",
    resources:
        runtime=7,
    log:
        logs + "/midway_metrics.{switch}",
    benchmark:
        bench + "/midway_metrics.{switch}",
    conda:
        "happler"
    shell:
        "workflow/scripts/midway_summary_metrics.py {params.use_flex_axes}-o {output.png} {params.metrics} &>{log}"


rule finemap_metrics:
    """summarize metrics from many executions of finemap_beta"""
    input:
        metrics = expand(
            rules.finemap_beta.output.metrics,
            beta=config["mode_attrs"]["beta"],
            sampsize=config["sample_size"],
            allow_missing=True
        ),
    params:
        metrics = lambda wildcards: expand(rules.finemap_beta.output.metrics, switch=wildcards.switch, allow_missing=True),
        use_flex_axes = lambda wildcards: "--use-flex-axes-limits " if wildcards.switch == "interact" else "",
    output:
        png=out + "/finemap_summary_metrics.{switch}.pdf",
    resources:
        runtime=7,
    log:
        logs + "/finemap_metrics.{switch}",
    benchmark:
        bench + "/finemap_metrics.{switch}",
    conda:
        "happler"
    shell:
        "workflow/scripts/midway_summary_metrics.py {params.use_flex_axes}-o {output.png} {params.metrics} &>{log}"
