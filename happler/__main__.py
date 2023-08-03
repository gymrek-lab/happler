#!/usr/bin/env python

import sys
import click
from pathlib import Path
from typing import TextIO
from typing import Union, Tuple

from . import tree
from haptools import data
from haptools import logging


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
def main():
    """
    happler: A haplotype-based fine-mapping method

    Test for associations between a trait and haplotypes (ie sets of correlated SNPs) rather than individual SNPs
    """
    pass


@main.command(context_settings=CONTEXT_SETTINGS)
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
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
    "-s",
    "--sample",
    "samples",
    type=str,
    multiple=True,
    show_default="all samples",
    help=(
        "A list of the samples to subset from the genotypes file (ex: '-s sample1 -s"
        " sample2')"
    ),
)
@click.option(
    "-S",
    "--samples-file",
    type=click.File("r"),
    show_default="all samples",
    help=(
        "A single column txt file containing a list of the samples (one per line) to"
        " subset from the genotypes file"
    ),
)
@click.option(
    "--discard-multiallelic",
    is_flag=True,
    show_default="do not discard multi-allelic variants",
    help="Whether to discard multi-allelic variants or just complain about them.",
)
@click.option(
    "--discard-missing",
    is_flag=True,
    show_default=True,
    default=False,
    help="Ignore any samples that are missing genotypes for the required variants",
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help="If using a PGEN file, read genotypes in chunks of X variants; reduces memory",
)
@click.option(
    "--maf",
    type=float,
    default=None,
    show_default="no filtering",
    help="Ignore variants with a MAF below this threshold",
)
@click.option(
    "--phased",
    is_flag=True,
    show_default=True,
    default=False,
    hidden=True,
    help="Do not check that variants are phased. Saves time and memory.",
)
@click.option(
    "-t",
    "--threshold",
    type=float,
    default=0.05,
    show_default=True,
    hidden=True,
    help="The alpha threshold used to determine when to terminate tree building",
)
@click.option(
    "--ld-prune-thresh",
    type=float,
    default=0.50,
    show_default=True,
    help="The LD threshold used to prune leaf nodes based on LD with their siblings",
)
@click.option(
    "--show-tree",
    is_flag=True,
    show_default=True,
    default=False,
    help="Output a tree in addition to the regular output.",
)
@click.option(
    "--covars",
    type=click.Path(exists=True, path_type=Path),
    show_default="no covariates",
    help="Any covariates to include in the model (as a PLINK2 .covar) file",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A .hap file describing the extracted haplotypes",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def run(
    genotypes: Path,
    phenotypes: Path,
    region: str = None,
    samples: Tuple[str] = tuple(),
    samples_file: Path = None,
    discard_multiallelic: bool = False,
    discard_missing: bool = False,
    chunk_size: int = None,
    maf: float = None,
    phased: bool = False,
    threshold: float = 0.05,
    ld_prune_thresh: float = None,
    show_tree: bool = False,
    covars: Path = None,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "INFO",
):
    """
    Use the tool to find trait-associated haplotypes

    GENOTYPES must be formatted as VCFs and

    PHENOTYPES must be a tab-separated file containing two columns: sample ID and
    phenotype value

    Ex: happler run tests/data/simple.vcf tests/data/simple.tsv > simple.hap
    """
    log = logging.getLogger(name="happler", level=verbosity)
    # handle samples
    if samples and samples_file:
        raise click.UsageError(
            "You may only use one of --sample or --samples-file but not both."
        )
    if samples_file:
        with samples_file as samps_file:
            samples = samps_file.read().splitlines()
    elif samples:
        # needs to be converted from tuple to list
        samples = list(samples)
    else:
        samples = None
    # load data
    log.info("Loading genotypes")
    if genotypes.suffix == ".pgen":
        gt = data.GenotypesPLINK(fname=genotypes, log=log, chunk_size=chunk_size)
    else:
        gt = data.GenotypesVCF(fname=genotypes, log=log)
    gt._prephased = phased
    gt.read(region=region, samples=samples)
    num_variants, num_samples = len(gt.variants), len(gt.samples)
    gt.check_missing(discard_also=discard_missing)
    removed = num_samples - len(gt.samples)
    if removed:
        log.info(f"Ignoring {removed} samples that are missing variants")
    gt.check_biallelic(discard_also=discard_multiallelic)
    removed = num_variants - len(gt.variants)
    if removed:
        log.info(f"Ignoring {removed} multiallelic variants")
        num_variants = len(gt.variants)
    gt.check_maf(threshold=maf, discard_also=True)
    removed = num_variants - len(gt.variants)
    if maf is not None:
        log.info(f"Ignoring {removed} variants with MAF < {maf}")
    gt.check_phase()
    log.info("There are {} samples and {} variants".format(*gt.data.shape))
    gt_samples_set = set(gt.samples)
    log.info("Loading phenotypes")
    ph = data.Phenotypes(fname=phenotypes, log=log)
    ph.read(samples=gt_samples_set)
    ph.standardize()
    ph.subset(samples=gt.samples, inplace=True)
    if len(ph.samples) < len(gt.samples):
        diff = gt_samples_set - set(ph.samples)
        log.error(
            f"The phenotypes file is missing {len(diff)} samples. Here are the first "
            f"few: {list(diff)[:5]}"
        )
    if len(ph.names) > 1:
        log.warning("Ignoring all but the first trait in the phenotypes file")
        ph.names = ph.names[:1]
        ph.data = ph.data[:, :1]
    log.info("Running tree builder")
    if covars:
        cv = data.Covariates(fname=covars, log=log)
        cv.read(samples=gt_samples_set)
        cv.subset(samples=gt.samples, inplace=True)
        if len(cv.samples) < len(gt.samples):
            diff = gt_samples_set - set(ph.samples)
            log.error(
                f"The covariates file is missing {len(diff)} samples. Here are the"
                f" first few: {list(diff)[:5]}"
            )
        test_method = tree.assoc_test.AssocTestSimpleCovariates(covars=cv.data)
    else:
        if ld_prune_thresh is None:
            test_method = tree.assoc_test.AssocTestSimpleSM()
        else:
            test_method = tree.assoc_test.AssocTestSimpleSMTScore()
    terminator = tree.terminator.TTestTerminator(thresh=threshold, log=log)
    hap_tree = tree.TreeBuilder(
        gt, ph, method=test_method, terminator=terminator, ld_prune_thresh=ld_prune_thresh, log=log
    ).run()
    log.info("Outputting haplotypes")
    tree.Haplotypes.from_tree(fname=output, tree=hap_tree, gts=gt, log=log).write()
    if show_tree:
        if output.suffix == ".gz":
            dot_output = output.with_suffix("").with_suffix(".dot")
        else:
            dot_output = output.with_suffix(".dot")
        log.info("Writing tree to dot file")
        with open(dot_output, "w") as dot_file:
            dot_file.write(hap_tree.dot())


if __name__ == "__main__":
    # run the CLI if someone tries 'python -m happler' on the command line
    main(prog_name="happler")
