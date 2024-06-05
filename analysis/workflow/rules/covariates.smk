from pathlib import Path


out = config["out"] + "/covars"
logs = out + "/logs"
bench = out + "/bench"


def check_config(value, default=False, place=config, as_set=False):
    """return true if config value exists and is true"""
    value = place[value] if (value in place and place[value]) else default
    return (set(value) if isinstance(value, list) else {value}) if as_set else value

rule peer:
    """compute PEER factors"""
    input:
        pheno_matrix = config["pheno_matrix"],
    output:
        covar=out+"/peer_factors.covar",
    resources:
        runtime=30,
        queue="hotel",
    log:
        logs + "/peer",
    benchmark:
        bench + "/peer",
    threads: 36
    conda:
        "../envs/peer.yml"
    shell:
        "workflow/scripts/peer_residuals.R {input} {output} &>{log}"


rule merge:
    """combine covariates into a merged .covar file"""
    input:
        covar = config["covar"],
        peers = rules.peer.output.covar,
    output:
        covar=out+"/covars.covar"
    resources:
        runtime=3,
    log:
        logs + "/covariates",
    benchmark:
        bench + "/covariates",
    conda:
        "../envs/default.yml"
    shell:
        "join -t $'\\t' -j1 "
        "<(sort -k1,1 {input.covar}) <(sort -k1,1 {input.peers}) > {output.covar}"
