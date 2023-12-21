# scripts
This directory contains various scripts used by the pipeline. However, you can use most of these scripts on their own, too. Some may even be helpful in day-to-day use.

All python scripts implement the --help argument. For bash and R scripts, you can run `head <script>` to read about their usage.

## [choose_different_ld.py](choose_different_ld.py)
Choose haplotype alleles such that their SNPs have a range of strengths of LD with the causal haplotype.

## [create_hap_file.sh](create_hap_file.sh)
Create a [.hap file](https://haptools.readthedocs.io/en/stable/formats/haplotypes.html) for use by haptools when simulating causal variables.

## [finemap_1p4.R](finemap_1p4.R)
Generate proper input files to run FINEMAP. Then run it.

## [heatmap_alleles.py](heatmap_alleles.py)
Create a heatmap of the phenotypes and alleles in a haplotype from a `.hap` file.

## [keep_samps.bash](keep_samps.bash)
Determine the intersection of the set of samples between multiple input reference panels

## [ld_heatmap.py](ld_heatmap.py)
Creates a heatmap of the LD pattern among the variants in a genotype matrix.

## [manhattan.py](manhattan.py)
Create a manhattan plot to visualize the results of a GWAS.

## [merge_plink.py](merge_plink.py)
Merge variants from two PGEN files that have the same set of samples.

## [parameter_plot.py](parameter_plot.py)
Plot the LD between a causal haplotype and a set of observed haplotypes across a range of hyper-parameters.

## [plot_phenotypes.py](plot_phenotypes.py)
Performs a linear regression and outputs a plot for phenotypes from `generate_phenotypes.py`.

## [run_FINEMAP.R](run_FINEMAP.R)
Run FINEMAP with a set of genotypes and phenotypes. Compute summary statistics, then call the `finemap_1p4.R` script.

## [run_SuSiE.R](run_SuSiE.R)
Run SuSiE with a set of genotypes and phenotypes. Follow the same interface as `run_FINEMAP.R`.

## [snakemake_io.py](snakemake_io.py)
Methods from the `snakemake.io` module which are imported by some of the other scripts in this directory.

## [summarize_results.R](summarize_results.R)
This R script summarizes the output of several executions of FINEMAP and SuSiE by producing several useful visualizations.

## [susie_metrics.R](susie_metrics.R)
Output specific metrics about the success of a SuSiE run.

## [utils.R](utils.R)
Define some useful methods for reading various file types into R. This script is sourced by other R scripts but never executed directly.
