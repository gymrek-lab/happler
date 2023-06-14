from pathlib import Path


out = config["out"] + "/happler"
logs = out + "/logs"
bench = out + "/bench"


rule run_happler:
    """ execute happler! """
    input:
        gts=config["snp_panel"],
        pts=config["pheno"],
    params:
        thresh=0.05,
    output:
        hap=out + "/happler/happler.hap",
        dot=out + "/happler/happler.dot",
    resources:
        runtime="0:15:00",
        queue="hotel",
    threads: 6
    log:
        logs + "/run_happler",
    benchmark:
        bench + "/run_happler",
    conda:
        "happler"
    shell:
        "happler run -o {output.hap} --verbosity DEBUG --discard-multiallelic"
        " -t {params.thresh} --show-tree {input.gts} {input.pts} &>{log}"


rule hap_tree:
    """ visualize the haplotype tree as a png file """
    input:
        dot=rules.run_happler.output.dot,
    params:
        file_ext = lambda wildcards, output: Path(output.png).suffix[1:],
    output:
        png=out + "/happler/happler.png",
    log:
        logs + "/hap_tree",
    benchmark:
        bench + "/hap_tree",
    conda:
        "../envs/default.yml"
    shell:
        "dot -T{params.file_ext} {input.dot} -o {output.png} &>{log}"


rule transform_happler:
    input:
        hap=rules.run_happler.output.hap,
        pgen=config["snp_panel"],
        pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        psam=Path(config["snp_panel"]).with_suffix(".psam"),
    params:
        hap_id="H1",
    output:
        pgen=temp(out + "/happler/happler.pgen"),
        pvar=temp(out + "/happler/happler.pvar"),
        psam=temp(out + "/happler/happler.psam"),
    resources:
        runtime="0:04:00"
    log:
        logs + "/transform_happler",
    benchmark:
        bench + "/transform_happler",
    conda:
        "happler"
    shell:
        "haptools transform -o {output.pgen} --id {params.hap_id} {input.pgen} "
        "{input.hap} &>{log}"


rule merge_happler:
    input:
        gts=config["snp_panel"],
        gts_pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        gts_psam=Path(config["snp_panel"]).with_suffix(".psam"),
        hps=rules.transform_happler.output.pgen,
        hps_pvar=rules.transform_happler.output.pvar,
        hps_psam=rules.transform_happler.output.psam,
    output:
        pgen=out + "/happler/merged_happler.pgen",
        pvar=out + "/happler/merged_happler.pvar",
        psam=out + "/happler/merged_happler.psam",
    resources:
        runtime="0:04:00"
    log:
        logs + "/merge_happler",
    benchmark:
        bench + "/merge_happler",
    conda:
        "happler"
    shell:
        "workflow/scripts/merge_plink.py {input.gts} {input.hps} {output.pgen} &> {log}"
