# scripts
This directory contains various scripts used by the pipeline. However, you can use most of these scripts on their own, too. Some may even be helpful in day-to-day use.

All python scripts implement the --help argument. For bash, awk, and R scripts, you can run `head <script>` to read about their usage.

## [check_resource_usage.py](check_resource_usage.py)
Using the log of a completed workflow run, compare existing memory and runtime allocations against benchmarked usage.

## [choose_different_ld.py](choose_different_ld.py)
Choose haplotype alleles such that their SNPs have a range of strengths of LD with the causal haplotype.

## [compute_pgen_ld.py](compute_pgen_ld.py)
Quickly compute the LD between a single variant and all other variants in a PGEN file.

## [conditional_regression_plots.py](conditional_regression_plots.py)
Make manhattan plots for a region after conditioning on 1) haplotypes and 2) their alleles together and 3) separately. Also, 4) show the original manhattan plot without any conditioning.

## [create_hap_file.sh](create_hap_file.sh)
Create a [.hap file](https://haptools.readthedocs.io/en/stable/formats/haplotypes.html) for use by haptools when simulating causal variables.

## [extract_pips.R](extract_pips.R)
Extract PIPs from the output of run_SuSiE.R

## [finemap_1p4.R](finemap_1p4.R)
Generate proper input files to run FINEMAP. Then run it.

## [heatmap_alleles.py](heatmap_alleles.py)
Create a heatmap of the phenotypes and alleles in a haplotype from a `.hap` file.

## [igv.bat](igv.bat)
Batch scripts for automated IGV plotting

## [keep_samps.bash](keep_samps.bash)
Determine the intersection of the set of samples between multiple input reference panels

## [ld_heatmap.py](ld_heatmap.py)
Creates a heatmap of the LD pattern among the variants in a genotype matrix.

## [manhattan.py](manhattan.py)
Create a manhattan plot to visualize the results of a GWAS.

## [merge_plink.py](merge_plink.py)
Merge variants from two PGEN files that have the same set of samples.

## [midway_manhattan_summary.py](midway_manhattan_summary.py)
Summarize the results of multiple midway manhattan plots via an ROC curve

## [midway_manhattan.bash](midway_manhattan.bash)
Create a "midway" manhattan plot depicting the p-values on a branch of a tree mid-way through happler's tree-building process.

## [parameter_plot.py](parameter_plot.py)
Plot the LD between a causal haplotype and a set of observed haplotypes across a range of hyper-parameters.

## [peer_residuals.R](peer_residuals.R)
Create PEER factors for gene expression phenotypes

## [plot_gts.py](plot_gts.py)
Plot phenotypes against genotypes for a haplotype.

## [plot_phenotypes.py](plot_phenotypes.py)
Performs a linear regression and outputs a plot for phenotypes from `generate_phenotypes.py`.

## [residuals.py](residuals.py)
Regresses out covariates from a set of phenotypes.

## [run_FINEMAP.R](run_FINEMAP.R)
Run FINEMAP with a set of genotypes and phenotypes. Compute summary statistics, then call the `finemap_1p4.R` script.

## [run_SuSiE.R](run_SuSiE.R)
Run SuSiE with a set of genotypes and phenotypes. Follow the same interface as `run_FINEMAP.R`.

## [snakemake_io.py](snakemake_io.py)
Methods from the `snakemake.io` module which are imported by some of the other scripts in this directory.

## [split_pheno.awk](split_pheno.awk)
Splits phenotypes from a .pheno file into multiple .pheno files.

## [summarize_results.R](summarize_results.R)
This R script summarizes the output of several executions of FINEMAP and SuSiE by producing several useful visualizations.

## [susie_metrics.R](susie_metrics.R)
Output specific metrics about the success of a SuSiE run.

## [utils.R](utils.R)
Define some useful methods for reading various file types into R. This script is sourced by other R scripts but never executed directly.

## [variance_explained_plot.py](variance_explained_plot.py)
Plot variance explained by haplotypes vs SNPs among all loci in a dataset.
