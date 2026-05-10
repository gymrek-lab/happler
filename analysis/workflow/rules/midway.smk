from pathlib import Path

out = config["out"]
logs = out + "/logs"
bench = out + "/bench"


# mapping from switch wildcard to integer for bash script
tswitch = {
    "interact": 1,
    "tscore" : 2,
    "covariance": 3,
    "bic": 4,
    "interact-bic": 5,
    "extension-bic": 6,
    "extension-tscore": 7,
}

wildcard_constraints:
    switch="("+"|".join(tswitch.keys())+")"

# if the pvar size is larger than 100 MB, just use the default memory instead
rsrc_func = lambda x: max if 100 > Path(x).with_suffix(".pvar").stat().st_size/1000/1000 else min


rule manhattan:
    """ run happler midway-through via a plink2 GWAS """
    input:
        gts=config["snp_panel"],
        gts_pvar=Path(config["snp_panel"]).with_suffix(".pvar"),
        gts_psam=Path(config["snp_panel"]).with_suffix(".psam"),
        pts=config["pheno"],
        hap=config["hap_file"],
    params:
        out_prefix = lambda wildcards, output: str(output.dir) + "/out",
        maf=config["min_maf"],
        target = "H0",
        tswitch=lambda wildcards: tswitch[wildcards.switch],
        just_target_snp=1,
        sim_mode=lambda wildcards: wildcards.sim_mode,
    output:
        dir=directory(out + "/{switch}"),
        linear=out + "/{switch}/out.linear",
        transform_pgen=temp(out + "/{switch}/out.pgen"),
        transform_pvar=temp(out + "/{switch}/out.pvar"),
        transform_psam=temp(out + "/{switch}/out.psam"),
        hap=out + "/{switch}/out.hap",
    resources:
        runtime=20,
        mem_mb=lambda wildcards, input: (
            rsrc_func(input.gts)(4000, Path(input.gts).with_suffix(".pvar").stat().st_size/1000 * 3.1342470426950246 + 1000)
        ) if str(wildcards.switch).startswith("extension") else 2000,
    log:
        logs + "/{switch}/manhattan",
    benchmark:
        bench + "/{switch}/manhattan",
    conda:
        "happler"
    shell:
        "mkdir -p {output.dir} && "
        "new_hap=\"$(workflow/scripts/flip_hap_alleles.py {input.gts} {input.pts} {input.hap} 2> {log})\" && "
        "{{ [ -z \"$new_hap\" ] || [ \"{params.sim_mode}\" != \"hap\" ] && cp {input.hap} {output.hap} || (echo \"$new_hap\" > {output.hap}); }} 2>>{log} && "
        "rsid=\"$(grep -E '^V' {output.hap} | cut -f5 | tail -n1)\" && "
        "workflow/scripts/midway_manhattan.bash {input.gts} {input.pts} {output.hap} "
        "{params.out_prefix} {params.target} \"$rsid\" {params.maf} {params.tswitch} "
        "{params.just_target_snp} &>> {log}"
