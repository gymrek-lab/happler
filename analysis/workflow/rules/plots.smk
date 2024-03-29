from pathlib import Path

out = config["out"]
logs = out + "/logs"
bench = out + "/bench"

mode = config["mode"]
out += "/" + mode
logs += "/" + mode
bench += "/" + mode

def agg_ld_range_obs(wildcards):
    """ return a list of hap files from happler """
    if "ld_range_checkpoint" in config:
        checkpoint_output = config["ld_range_checkpoint"].get(**wildcards).output.hap
        return expand(
            config["happler_hap"],
            ld=glob_wildcards(Path(checkpoint_output) / "ld_{ld}/haplotype.hap").ld,
            beta=config["mode_attrs"]["beta"],
        )
    else:
        return expand(config["happler_hap"], beta=config["mode_attrs"]["beta"])

def agg_ld_range_causal(wildcards):
    """ return a list of hap files from the LD range checkpoint """
    if "ld_range_checkpoint" in config:
        checkpoint_output = config["ld_range_checkpoint"].get(**wildcards).output.hap
        return expand(
            str(checkpoint_output),
            ld=glob_wildcards(Path(checkpoint_output) / "ld_{ld}/haplotype.hap").ld,
            beta=config["mode_attrs"]["beta"],
        )
    else:
        causal_hap = config["causal_hap"]
        return expand(causal_hap, beta=config["mode_attrs"]["beta"])


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
        runtime="0:10:00",
    log:
        logs + "/plot_params",
    benchmark:
        bench + "/plot_params",
    conda:
        "happler"
    shell:
        "workflow/scripts/parameter_plot.py -o {output.png} "
        "{input.gts} {params.observed_haps} {params.causal_hap} &> {log}"
