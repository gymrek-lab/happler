from pathlib import Path


out = config["out"]
logs = out + "/logs"
bench = out + "/bench"


rule run:
    """ execute happler! """
    input:
        gts=config["snp_panel"],
        pts=config["pheno"],
    params:
        thresh=0.05,
    output:
        hap=out + "/happler.hap",
        dot=out + "/happler.dot",
    resources:
        runtime="0:30:00",
        queue="hotel",
    threads: 6
    log:
        logs + "/run",
    benchmark:
        bench + "/run",
    conda:
        "happler"
    shell:
        "happler run -o {output.hap} --verbosity DEBUG "
        "--discard-multiallelic "
        "-t {params.thresh} --show-tree {input.gts} {input.pts} &>{log}"


rule tree:
    """ visualize the haplotype tree as a png file """
    input:
        dot=rules.run.output.dot,
    params:
        file_ext = lambda wildcards, output: Path(output.png).suffix[1:],
    output:
        png=out + "/happler.png",
    log:
        logs + "/tree",
    benchmark:
        bench + "/tree",
    conda:
        "../envs/default.yml"
    shell:
        "dot -T{params.file_ext} {input.dot} -o {output.png} &>{log}"


rule transform:
    input:
        hap=rules.run.output.hap,
        pgen=config["snp_panel"],
        pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        psam=Path(config["snp_panel"]).with_suffix(".psam"),
    params:
        hap_id="H1",
    output:
        pgen=temp(out + "/happler.pgen"),
        pvar=temp(out + "/happler.pvar"),
        psam=temp(out + "/happler.psam"),
    resources:
        runtime="0:04:00"
    log:
        logs + "/transform",
    benchmark:
        bench + "/transform",
    conda:
        "happler"
    shell:
        "haptools transform -o {output.pgen} --id {params.hap_id} {input.pgen} "
        "{input.hap} &>{log}"


rule merge:
    input:
        gts=config["snp_panel"],
        gts_pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        gts_psam=Path(config["snp_panel"]).with_suffix(".psam"),
        hps=rules.transform.output.pgen,
        hps_pvar=rules.transform.output.pvar,
        hps_psam=rules.transform.output.psam,
    output:
        pgen=out + "/merged.pgen",
        pvar=out + "/merged.pvar",
        psam=out + "/merged.psam",
    resources:
        runtime="0:04:00"
    log:
        logs + "/merge",
    benchmark:
        bench + "/merge",
    conda:
        "happler"
    shell:
        "workflow/scripts/merge_plink.py {input.gts} {input.hps} {output.pgen} &> {log}"


rule finemapper:
    """ execute SuSiE using the haplotypes from happler """
    input:
        gt=rules.merge.output.pgen,
        phen=config["pheno"],
    params:
        outdir=lambda wildcards, output: Path(output.susie).parent,
        exclude_causal="NULL",
    output:
        susie=out + "/susie.rds",
    resources:
        runtime="1:15:00",
        queue="hotel",
    log:
        logs + "/finemapper",
    benchmark:
        bench + "/finemapper",
    conda:
        "../envs/susie.yml"
    shell:
        "workflow/scripts/run_SuSiE.R {input} {params} &>{log}"


rule results:
    """
        create plots to summarize the results of the simulations when tested
        on happler
    """
    input:
        gt=rules.finemapper.input.gt,
        susie=rules.finemapper.output.susie,
        happler_hap=rules.run.output.hap,
    params:
        outdir=lambda wildcards, output: Path(output.susie_pdf).parent,
        exclude_causal=lambda wildcards: 0,
        causal_hap=config["hap_file"],
        causal_gt=config["causal_gt"],
    output:
        susie_pdf = out + "/susie.pdf",
        # susie_eff_pdf=temp(out + "/susie_eff.pdf"),
    log:
        logs + "/results",
    conda:
        "../envs/susie.yml"
    script:
        "../scripts/summarize_results.R"


rule gwas:
    """run a GWAS"""
    input:
        pgen=rules.merge.output.pgen,
        pvar=rules.merge.output.pvar,
        psam=rules.merge.output.psam,
        pts=config["pheno"],
    params:
        in_prefix = lambda w, input: Path(input.pgen).with_suffix(""),
        out_prefix = lambda w, output: Path(output.log).with_suffix(""),
    output:
        log = temp(out + "/hap.log"),
        linear = out + "/hap.hap.glm.linear",
    resources:
        runtime="0:10:00",
    log:
        logs + "/gwas",
    benchmark:
        bench + "/gwas",
    threads: 1
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --linear allow-no-covars --variance-standardize "
        "--pheno iid-only {input.pts} --pfile {params.in_prefix} "
        "--out {params.out_prefix} --threads {threads} &>{log}"


rule manhattan:
    input:
        linear=rules.gwas.output.linear,
    params:
        linear = lambda wildcards, input: f"-l "+input.linear,
        red_ids = lambda wildcards: [
            f"-i {i.split(':')[0]}" for i in config["snps"]
        ],
        orange_ids = lambda wildcards: "-b hap -b H1",
    output:
        png = out + "/manhattan.pdf",
    resources:
        runtime="0:05:00"
    log:
        logs + "/manhattan",
    benchmark:
        bench + "/manhattan",
    conda:
        "happler"
    shell:
        "workflow/scripts/manhattan.py -o {output.png} {params.linear} "
        "{params.red_ids} {params.orange_ids} &>{log}"
