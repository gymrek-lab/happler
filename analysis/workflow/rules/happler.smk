from pathlib import Path


out = config["out"]
logs = out + "/logs"
bench = out + "/bench"


# whether to include the happler haplotype ("in") or no haplotypes ("ex")
exclude_obs = {"in": 0, "ex": 1}
# or, if config["random"] is not None, this denotes
# whether to include a random haplotype ("in") or the causal haplotype ("ex")

def check_config(value, default=False, place=config, as_set=False):
    """return true if config value exists and is true"""
    value = place[value] if (value in place and place[value]) else default
    return (set(value) if isinstance(value, list) else {value}) if as_set else value

def parse_locus(locus):
    """parse locus into chrom, start, end"""
    chrom = locus.split("_")[0]
    end = locus.split("-")[1]
    start = locus.split("_")[1].split("-")[0]
    return chrom, start, end


rule run:
    """ execute happler! """
    input:
        gts=config["snp_panel"],
        pts=config["pheno"],
        covar=config["covar"],
    params:
        thresh=lambda wildcards: 0.05 if "alpha" not in wildcards else wildcards.alpha,
        region=lambda wildcards: wildcards.locus.replace("_", ":"),
        covar=lambda wildcards, input: ("--covar " + input["covar"] + " ") if check_config("covar") else "",
        maf = check_config("min_maf", 0),
    output:
        hap=out + "/happler.hap",
        gz=out + "/happler.hap.gz",
        idx=out + "/happler.hap.gz.tbi",
        dot=out + "/happler.dot",
    resources:
        runtime=10,
        # slurm_partition="hotel",
        # slurm_extra="--qos=hotel",
        # mem_mb=lambda wildcards, threads: threads*4.57,
    threads: 6
    log:
        logs + "/run",
    benchmark:
        bench + "/run",
    conda:
        "happler"
    shell:
        "happler run -o {output.hap} --verbosity DEBUG --maf {params.maf} "
        "--discard-multiallelic --region {params.region} {params.covar}"
        "-t {params.thresh} --show-tree {input.gts} {input.pts} &>{log} && "
        "haptools index -o {output.gz} {output.hap} &>>{log}"


rule tree:
    """ visualize the haplotype tree as a png file """
    input:
        dot=rules.run.output.dot,
    params:
        file_ext = lambda wildcards, output: Path(output.png).suffix[1:],
    output:
        png=out + "/happler.png",
    resources:
        runtime=4,
    log:
        logs + "/tree",
    benchmark:
        bench + "/tree",
    conda:
        "../envs/default.yml"
    shell:
        "dot -T{params.file_ext} {input.dot} -o {output.png} &>{log}"


rule cond_linreg:
    """plot conditional regressions for a haplotype"""
    input:
        pgen=config["snp_panel"],
        pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        psam=Path(config["snp_panel"]).with_suffix(".psam"),
        hap=rules.run.output.hap,
        pts=config["pheno"],
    params:
        hap_id = "H0",
        maf = check_config("min_maf", 0),
        region=lambda wildcards: wildcards.locus.replace("_", ":"),
    output:
        png=out + "/cond_linreg.pdf",
    resources:
        runtime=15,
    log:
        logs + "/cond_linreg",
    benchmark:
        bench + "/cond_linreg",
    conda:
        "happler"
    shell:
        "workflow/scripts/conditional_regression_plots.py --verbosity DEBUG "
        "--show-original --region {params.region} -o {output.png} -i {params.hap_id} "
        "--maf {params.maf} {input.pgen} {input.pts} {input.hap} &>{log}"


rule heatmap:
    """look at the LD pattern of the haplotype"""
    # TODO: also include causal hap if one exists
    input:
        pgen=config["snp_panel"],
        pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        psam=Path(config["snp_panel"]).with_suffix(".psam"),
        hap=rules.run.output.hap,
        pts=config["pheno"],
    params:
        region=lambda wildcards: wildcards.locus.replace("_", ":"),
    output:
        png=out + "/heatmap.pdf",
    resources:
        runtime=4,
    log:
        logs + "/heatmap",
    benchmark:
        bench + "/heatmap",
    conda:
        "happler"
    shell:
        "workflow/scripts/heatmap_alleles.py --verbosity DEBUG "
        "--region {params.region} -o {output.png} "
        "{input.pgen} {input.hap} {input.pts} &>{log}"


rule igv:
    """look at the haplotype in IGV"""
    input:
        hap=rules.run.output.hap,
        prefs = "data/igv.prefs.properties", # todo move this into config
        hprc_urls = "data/hprc.bam.list",
    params:
        chrom=lambda wildcards: parse_locus(wildcards.locus)[0],
        outdir=lambda wildcards, output: Path(output.png).parent,
        outfile=lambda wildcards, output: Path(output.png).name,
        hap_id="H0", #you can change this to "" to use all haps
        width=10000,
    output:
        bed=temp(out + "/happler_snps.bed"),
        png=out + "/igv.png",
    resources:
        runtime=10,
    log:
        logs + "/igv",
    benchmark:
        bench + "/igv",
    conda:
        "../envs/default.yml"
    shell:
        "region={params.chrom}:$("
        "grep -Ev '^#' {input.hap} | awk -F $'\\t' "
        "'BEGIN {{ min = \"inf\"; max = \"-inf\" }} "
        "'\"$([ \"{params.hap_id}\" != \"\" ] && echo '$2 == \"{params.hap_id}\" ')\"'"
        "{{ if ($3 < min) min = $3; if ($4 > max) max = $4 }} END {{ print (min-{params.width})\"-\"(max+{params.width}) }}'"
        ") && grep -E '^V' {input.hap} | "
        "awk -F $'\\t' -v 'OFS=\\t' "
        "\"$([ \"{params.hap_id}\" != \"\" ] && echo '$2 == \"{params.hap_id}\" ')\"'"
        "{{print {params.chrom}, $3, $4, $2\"_\"$5;}}' > {output.bed} && "
        "igv -o {input.prefs} -b <("
        "sed 's+mySnapshotDirectory+{params.outdir}+;s/REGION/'$region'/;s+BED+{output.bed}+' workflow/scripts/igv.bat && "
        "cut -f1,2 --output-delimiter=: {output.bed} | sort -u | sed 's/^/SORT BASE /' && echo 'snapshot {params.outfile}' && echo exit"
        ") &>{log}"


rule transform:
    input:
        hap=rules.run.output.gz,
        pgen=config["snp_panel"],
        pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        psam=Path(config["snp_panel"]).with_suffix(".psam"),
    params:
        region=lambda wildcards: wildcards.locus.replace("_", ":"),
    output:
        pgen=temp(out + "/happler.pgen"),
        pvar=temp(out + "/happler.pvar"),
        psam=temp(out + "/happler.psam"),
    resources:
        runtime=4,
    log:
        logs + "/transform",
    benchmark:
        bench + "/transform",
    conda:
        "happler"
    shell:
        "haptools transform -o {output.pgen} --region {params.region} "
        "{input.pgen} {input.hap} &>{log}"


rule linreg:
    """plot a haplotype's genotypes against its phenotypes"""
    input:
        hap=rules.transform.output.pgen,
        hap_pvar=rules.transform.output.pvar,
        hap_psam=rules.transform.output.psam,
        pts=config["pheno"],
    params:
        region=lambda wildcards: wildcards.locus.replace("_", ":"),
    output:
        png=out + "/linreg.pdf",
    resources:
        runtime=4,
    log:
        logs + "/linreg",
    benchmark:
        bench + "/linreg",
    conda:
        "happler"
    shell:
        "workflow/scripts/plot_gts.py --verbosity DEBUG "
        "-o {output.png} --region {params.region} "
        "{input.hap} {input.pts} 2>{log}"


rule sv_ld:
    """compute LD between this haplotype and a bunch of SVs"""
    input:
        hap=rules.transform.output.pgen,
        hap_pvar=rules.transform.output.pvar,
        hap_psam=rules.transform.output.psam,
        sv=lambda wildcards: config["SVs"],
        sv_pvar=lambda wildcards: Path(config["SVs"]).with_suffix(".pvar"),
        sv_psam=lambda wildcards: Path(config["SVs"]).with_suffix(".psam"),
    params:
        start=lambda wildcards: max(0, int(parse_locus(wildcards.locus)[1])-1000000),
        end=lambda wildcards: int(parse_locus(wildcards.locus)[2])+1000000,
        chrom=lambda wildcards: parse_locus(wildcards.locus)[0],
        maf = check_config("min_maf", 0),
        hapid="H0",
    output:
        ld=out + "/happler_svs.ld",
    resources:
        runtime=4,
    log:
        logs + "/sv_ld",
    benchmark:
        bench + "/sv_ld",
    conda:
        "happler"
    shell:
        "workflow/scripts/compute_pgen_ld.py --verbosity DEBUG "
        "--region '{params.chrom}:{params.start}-{params.end}' --maf {params.maf} "
        "--hap-id {params.hapid} -o /dev/stdout {input.sv} {input.hap} 2>{log} | "
        "grep -Ev 'nan$' > {output} 2>>{log}"


def merge_hps_input(wildcards):
    if config["random"] is None:
        # include the hap that happler found
        return rules.transform.output
    else:
        if exclude_obs[wildcards.ex]:
            # exclude the random hap (and use the causal hap, instead)
            return config["causal_gt"]
        else:
            # include the random hap
            return config["random"]


rule merge:
    input:
        gts=config["snp_panel"],
        gts_pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        gts_psam=Path(config["snp_panel"]).with_suffix(".psam"),
        hps=lambda wildcards: merge_hps_input(wildcards).pgen,
        hps_pvar=lambda wildcards: merge_hps_input(wildcards).pvar,
        hps_psam=lambda wildcards: merge_hps_input(wildcards).psam,
    params:
        region=lambda wildcards: wildcards.locus.replace("_", ":"),
        maf = check_config("min_maf", 0),
    output:
        pgen=temp(out + "/{ex}clude/merged.pgen"),
        pvar=temp(out + "/{ex}clude/merged.pvar"),
        psam=temp(out + "/{ex}clude/merged.psam"),
    resources:
        runtime=10,
    log:
        logs + "/{ex}clude/merge",
    benchmark:
        bench + "/{ex}clude/merge",
    conda:
        "happler"
    shell:
        "workflow/scripts/merge_plink.py --region {params.region} --maf {params.maf} "
        "{input.gts} {input.hps} {output.pgen} &> {log}"


finemapper_input = lambda wildcards: rules.transform.input if (exclude_obs[wildcards.ex] and config["random"] is None) else rules.merge.output


rule finemapper:
    """ execute SuSiE using the haplotypes from happler """
    input:
        gt=lambda wildcards: finemapper_input(wildcards).pgen,
        gt_pvar=lambda wildcards: finemapper_input(wildcards).pvar,
        gt_psam=lambda wildcards: finemapper_input(wildcards).psam,
        phen=config["pheno"],
    params:
        outdir=lambda wildcards, output: Path(output.susie).parent,
        exclude_causal="NULL",
        region=lambda wildcards: wildcards.locus.replace("_", ":"),
    output:
        susie=out + "/{ex}clude/susie.rds",
    resources:
        runtime=75,
    threads: 6,
    log:
        logs + "/{ex}clude/finemapper",
    benchmark:
        bench + "/{ex}clude/finemapper",
    conda:
        "../envs/susie.yml"
    shell:
        "workflow/scripts/run_SuSiE.R {input.gt} {input.phen} {params} &>{log}"


rule pips:
    """ extract PIPs to a TSV file """
    input:
        rds=rules.finemapper.output.susie,
    output:
        tsv=out + "/{ex}clude/susie_pips.tsv",
    resources:
        runtime=5,
    threads: 3,
    log:
        logs + "/{ex}clude/pips",
    benchmark:
        bench + "/{ex}clude/pips",
    conda:
        "../envs/susie.yml"
    shell:
        "workflow/scripts/extract_pips.R {input} {output} &>{log}"


rule metrics:
    """ compute summary metrics from the output of the finemapper """
    input:
        finemap=expand(rules.finemapper.output.susie, ex="in", allow_missing=True),
        obs_hap=rules.run.output.hap,
        caus_hap=config["hap_file"],
    output:
        metrics=out + "/susie_metrics.tsv",
    resources:
        runtime=5,
    threads: 6,
    log:
        logs + "/metrics",
    benchmark:
        bench + "/metrics",
    conda:
        "../envs/susie.yml"
    shell:
        "workflow/scripts/susie_metrics.R {input} >{output} 2>{log}"


def results_happler_hap_input(wildcards):
    if config["random"] is None:
        if exclude_obs[wildcards.ex]:
            return []
        return rules.run.output.hap
    elif exclude_obs[wildcards.ex]:
        return expand(config["hap_file"], **wildcards)
    return expand(config["random_hap"], **wildcards)


def results_causal_hap_input(wildcards):
    if config["random"] is None:
        if exclude_obs[wildcards.ex]:
            return []
    return expand(config["hap_file"], **wildcards)


rule results:
    """
        create plots to summarize the results of the simulations when tested
        on happler
    """
    input:
        gt=rules.finemapper.input.gt,
        susie=rules.finemapper.output.susie,
        happler_hap=results_happler_hap_input,
        causal_gt=config["causal_gt"].pgen if "causal_gt" in config else [],
        causal_hap=results_causal_hap_input,
    params:
        outdir=lambda wildcards, output: Path(output.susie_pdf).parent,
        region=lambda wildcards: wildcards.locus.replace("_", ":"),
    output:
        susie_pdf = out + "/{ex}clude/susie.pdf",
        # susie_eff_pdf=temp(out + "/susie_eff.pdf"),
    resources:
        runtime=5,
    log:
        logs + "/{ex}clude/results",
    benchmark:
        bench + "/{ex}clude/results",
    conda:
        "../envs/susie.yml"
    script:
        "../scripts/summarize_results.R"


rule merge_SVs:
    input:
        gts=rules.merge.output.pgen,
        gts_pvar=rules.merge.output.pvar,
        gts_psam=rules.merge.output.psam,
        svs=lambda wildcards: config["SVs"],
        svs_pvar=lambda wildcards: Path(config["SVs"]).with_suffix(".pvar"),
        svs_psam=lambda wildcards: Path(config["SVs"]).with_suffix(".psam"),
    params:
        maf = check_config("min_maf", 0),
        region=lambda wildcards: wildcards.locus.replace("_", ":"),
    output:
        pgen=temp(out + "/{ex}clude/merged_SVs.pgen"),
        pvar=temp(out + "/{ex}clude/merged_SVs.pvar"),
        psam=temp(out + "/{ex}clude/merged_SVs.psam"),
    resources:
        runtime=4,
    log:
        logs + "/{ex}clude/merge_SVs",
    benchmark:
        bench + "/{ex}clude/merge_SVs",
    conda:
        "happler"
    shell:
        "workflow/scripts/merge_plink.py --region {params.region} --maf {params.maf} "
        "{input.gts} {input.svs} {output.pgen} &> {log}"


rule gwas:
    """run a GWAS"""
    input:
        pgen=(rules.merge_SVs if check_config("SVs") else rules.merge).output.pgen,
        pvar=(rules.merge_SVs if check_config("SVs") else rules.merge).output.pvar,
        psam=(rules.merge_SVs if check_config("SVs") else rules.merge).output.psam,
        pts=config["pheno"],
        covar=config["covar"],
    params:
        maf = check_config("min_maf", 0),
        in_prefix = lambda w, input: Path(input.pgen).with_suffix(""),
        out_prefix = lambda w, output: Path(output.log).with_suffix(""),
        covar = lambda wildcards, input: ("'no-x-sex' --covar 'iid-only' " + input["covar"]) if check_config("covar") else "allow-no-covars",
        start=lambda wildcards: parse_locus(wildcards.locus)[1],
        end=lambda wildcards: parse_locus(wildcards.locus)[2],
        chrom=lambda wildcards: parse_locus(wildcards.locus)[0],
    output:
        log = temp(out + "/{ex}clude/hap.log"),
        linear = out + "/{ex}clude/hap.hap.glm.linear",
    resources:
        runtime=10,
    log:
        logs + "/{ex}clude/gwas",
    benchmark:
        bench + "/{ex}clude/gwas",
    threads: 1
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --glm {params.covar} --variance-standardize --maf {params.maf} "
        "--pheno iid-only {input.pts} --pfile {params.in_prefix} "
        # "--from-bp {params.start} --to-bp {params.end} --chr {params.chrom} " # unnecessary b/c merge subsets by region already
        "--out {params.out_prefix} --threads {threads} &>{log}"


rule manhattan:
    input:
        linear=rules.gwas.output.linear,
    params:
        linear = lambda wildcards, input: f"-l "+input.linear,
        red_ids = lambda wildcards: [
            f"-i {i.split(':')[0]}" for i in check_config("snps", default=[])
        ] if not check_config("SVs") else "--red-chr-ids",
        orange_ids = lambda wildcards: "-b hap --orange-Hids",
        # TODO: allow specifying IDs from hap files, instead
    output:
        png = out + "/{ex}clude/manhattan.pdf",
    resources:
        runtime=10,
    log:
        logs + "/{ex}clude/manhattan",
    benchmark:
        bench + "/{ex}clude/manhattan",
    conda:
        "happler"
    shell:
        "workflow/scripts/manhattan.py -o {output.png} --no-label {params.linear} "
        "{params.red_ids} {params.orange_ids} &>{log}"
