from pathlib import Path

out = config["out"]
logs = out + "/logs"
bench = out + "/bench"


wildcard_constraints:
    switch="(interact|tscore|covariance|bic|interact-bic)"

# mapping from switch wildcard to integer for bash script
tswitch = {
    "interact": 1,
    "tscore" : 2,
    "covariance": 3,
    "bic": 4,
    "interact-bic": 5,
}


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
    output:
        dir=directory(out + "/{switch}"),
        linear=out + "/{switch}/out.linear",
        transform_pgen=temp(out + "/{switch}/out.pgen"),
        transform_pvar=temp(out + "/{switch}/out.pvar"),
        transform_psam=temp(out + "/{switch}/out.psam"),
    resources:
        runtime=10,
    log:
        logs + "/{switch}/manhattan",
    benchmark:
        bench + "/{switch}/manhattan",
    conda:
        "happler"
    shell:
        "rsid=\"$(grep -E '^V' {input.hap} | cut -f5 | tail -n1)\" && "
        "workflow/scripts/midway_manhattan.bash {input.gts} {input.pts} {input.hap} "
        "{params.out_prefix} {params.target} \"$rsid\" {params.maf} {params.tswitch} "
        "{params.just_target_snp} &> {log}"
