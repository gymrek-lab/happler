#!/usr/bin/env python
from pathlib import Path

# Allow us to edit fonts in Illustrator
import click
import matplotlib
import numpy as np
matplotlib.use('Agg')
from haptools import logging
import matplotlib.pyplot as plt
from haptools.data import GenotypesVCF, GenotypesPLINK, Haplotypes, Phenotypes


CAUSAL_COLOR_KEY = {
    "blue": "right SNP, right allele",
    "purple": "right SNP, wrong allele",
    "red": "observed wrong SNP",
    "black": "unobserved causal SNP",
}


def plot_hapmatrix(ax, hpmt, hap_id, snps, colors = None):
    """
    Adapted from this script by Shubham Saini
    https://github.com/shubhamsaini/pgc_analysis/blob/master/viz-snp-hap.py
    """
    # box_w =  1.0/hpmt.shape[1]
    # box_h = box_w
    # hap_height = hpmt.shape[0]*0.0025*4
    # legend_height = 0
    # Plot SNPs
    ax.imshow(hpmt, vmin=0, vmax=1, cmap=plt.cm.Greys, aspect="auto", interpolation="none")
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xticks(np.arange(0, len(snps), 1))
    ax.set_xticklabels(snps, fontsize=6.5)
    # map between list of indices and list of colors from CAUSAL_COLOR_KEY
    if colors is not None:
        colors = map(list(CAUSAL_COLOR_KEY.keys()).__getitem__, colors)
        for xtick, color in zip(ax.get_xticklabels(), colors):
            xtick.set_color(color)
    ax.set_title("All haplotypes" if hap_id == "ALL" else "Haplotype %s"%hap_id)


def get_cmap_string(palette, domain):
    """
    Uniquely map a list of strings to distinct colors

    Return a dictionary mapping each domain string to a color

    Adapted from https://stackoverflow.com/a/59034585/16815703
    """
    domain_unique = np.unique(domain)
    hash_table = {key: i_str for i_str, key in enumerate(domain_unique)}
    mpl_cmap = matplotlib.cm.get_cmap(palette, lut=len(domain_unique))

    def cmap_out(X, **kwargs):
        return mpl_cmap(hash_table[X], **kwargs)

    return dict(zip(domain, map(cmap_out, domain)))


def plot_hap_label_table(ax, hps, hps_vars):
    """
    Append a haplotype "label table" to the bottom of the haplotype matrix

    This is useful when multiple haplotypes are pictured in the table
    """
    hps_ids = list(hps.data.keys())
    # create a dictionary mapping haplotypes to unique colors
    hps_colors = get_cmap_string('viridis', hps_ids)
    # create a matrix of colors for the table
    snp_colors = np.array([
        [(hps_colors[hp.id] if snp in hp.varIDs else (1, 1, 1, 1)) for snp in hps_vars]
        for hp in hps.data.values()
    ], dtype='float32')
    # add a blank column for pheno
    snp_colors = np.append(
        snp_colors,
        np.ones((snp_colors.shape[0], 1, snp_colors.shape[2])),
        1
    )
    # now finally create the plot
    table = ax.table(
        cellColours=snp_colors,
        rowLabels=hps_ids,
        bbox=[0, -0.28, 1, 0.05*len(hps_ids)],
        # edges="horizontal",
        loc="bottom",
    )
    # remove table edges
    for c in table.properties()['celld'].values():
        c.set(linewidth=0)


@click.command()
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("phenotypes", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all genotypes",
    help="""
    The region from which to extract genotypes; ex: 'chr1:1234-34566' or 'chr7'\n
    For this to work, the VCF must be indexed and the seqname must match!""",
)
@click.option(
    "-i",
    "--hap-id",
    type=str,
    show_default="all of the haplotypes",
    help=(
        "A haplotype ID from the .hap file to plot"
        "(ex: '-i H1')."
    ),
)
@click.option(
    "-c",
    "--causal",
    type=click.Path(path_type=Path, exists=True),
    default=None,
    show_default=True,
    help=(
        "The path to a .hap file containing haplotypes simulated to be causal. "
        "These will be added to the plot and highlited differently."
    ),
)
@click.option(
    "--use-hap-alleles",
    is_flag=True,
    default=False,
    show_default=True,
    help="Color alleles from the haplotype as black instead of just the ALT alleles"
)
@click.option(
    "--label-haps/--no-label-haps",
    is_flag=True,
    default=True,
    show_default=True,
    help="Whether to add labels for each of the haplotypes. Doesn't work with --causal",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A PNG file containing the desired heatmap",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def main(
    genotypes: Path,
    haplotypes: Path,
    phenotypes: Path,
    region: str = None,
    hap_id: str = None,
    causal: Path = None,
    use_hap_alleles: bool = False,
    label_haps: bool = True,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "INFO"
):
    """
    Create a heatmap plot that visualizes a haplotype from a .hap file
    """
    log = logging.getLogger("heatmap_alleles", level=verbosity)

    gts = GenotypesVCF
    if genotypes.suffix == ".pgen":
        gts = GenotypesPLINK

    gts = gts(genotypes, log=log)

    hps = Haplotypes(haplotypes, log=log)
    hps.read(haplotypes=(set((hap_id,)) if hap_id is not None else None))

    # get the variants from all haplotypes
    hps_vars = {var.id: var.allele for hap in hps.data for var in hps.data[hap].variants}
    if hap_id is None:
        hap_id = "ALL"
        if use_hap_alleles:
            log.warning("A hap ID wasn't given. Setting --use-hap-alleles to False")
            use_hap_alleles = False

    snp_colors = None
    if causal is not None:
        # since --label-haps is not yet supported with --causal
        label_haps = False

        hps_causal = Haplotypes(causal, log=log)
        hps_causal.read()
        # if there is more than one haplotype
        if len(hps_causal.data) != 1:
            raise ValueError("The causal .hap file must contain only one haplotype")
        # extract causal alleles
        causal_hap_id = list(hps_causal.data.keys())[0]
        causal_hps_vars = {
            var.id: var.allele
            for var in hps_causal.data[causal_hap_id].variants
        }
        # mark whether observed alleles are correct/incorrect
        snp_colors = {
            var_id: (
                int(allele != causal_hps_vars[var_id])
                if var_id in causal_hps_vars else 2
            )
            for var_id, allele in hps_vars.items()
        }
        # add unobserved causal alleles to hps_vars and snp_colors
        for var_id in (causal_hps_vars.keys() - hps_vars.keys()):
            hps_vars[var_id] = causal_hps_vars[var_id]
            snp_colors[var_id] = 3

    pts = Phenotypes(phenotypes, log=log)
    pts.read()

    gts.read(variants=set(hps_vars.keys()), region=region, samples=pts.samples)
    gts.check_phase()
    gts.check_missing()
    gts.check_biallelic()
    gts.subset(variants=list(hps_vars.keys()), inplace=True)
    pts.subset(samples=gts.samples)
    num_samps, num_snps, _ = gts.data.shape

    if use_hap_alleles:
        for var_idx, var in enumerate(gts.variants):
            # is it the REF or ALT allele?
            allele = var["alleles"][1] == hps_vars[var["id"]]
            if snp_colors is not None and (snp_colors[var["id"]] == 1):
                allele = not allele
            # modify the genotype matrix to contain 1s only when the allele is present
            gts.data[:, var_idx] = gts.data[:, var_idx] == allele
            gts.variants[var_idx]["id"] += f":{int(allele)}"

    # resulting shape is num_samps * 2 by num_SNPs
    hpmt = gts.data.transpose((0, 2, 1)).reshape((num_samps*2, num_snps)).astype(np.bool_)
    # sort by each SNP from left to right
    samp_indices = np.lexsort(tuple(hpmt[:, i] for i in range(hpmt.shape[1]-1, -1, -1)))
    hpmt = hpmt[samp_indices]

    # also append the phenotypes
    pts = np.repeat(pts.data[:, 0], 2)[samp_indices]
    # standardize so that the max value is 1
    pts = (pts-pts.min())/(pts.max()-pts.min())
    hpmt = np.append(hpmt, pts[:, np.newaxis], axis=1)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_hapmatrix(
        ax, hpmt, hap_id, list(gts.variants["id"])+["pheno"], (
            snp_colors.values() if snp_colors is not None else None
        ),
    )
    if label_haps:
        plot_hap_label_table(ax, hps, list(gts.variants["id"]))
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig(output)


if __name__ == "__main__":
    main()
