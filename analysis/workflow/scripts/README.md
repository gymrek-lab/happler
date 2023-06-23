# scripts
This directory contains various scripts used by the pipeline. However, you can use most of these scripts on their own, too. Some may even be helpful in day-to-day use.

All python scripts implement the --help argument. For bash and R scripts, you can run `head <script>` to read about their usage.

## [choose_different_ld.py](choose_different_ld.py)
Choose haplotype alleles such that their SNPs have a range of strengths of LD with the causal haplotype.

## [create_hap_file.sh](create_hap_file.sh)
Create a [.hap file](https://haptools.readthedocs.io/en/stable/formats/haplotypes.html) for use by haptools when simulating causal variables.

## [finemap_1p4.R](finemap_1p4.R)
Generate proper input files to run FINEMAP. Then run it.

## [finemapping_methods.R](finemapping_methods.R)
Run FINEMAP and SuSiE. Duplicates code from `run_FINEMAP.R` and `run_SuSiE.R`.

## [generate_phenotypes.py](generate_phenotypes.py)
Uses pre-existing gentoype matrices to create a file of phenotypes.

## [keep_samps.bash](keep_samps.bash)
Determine the intersection of the set of samples between multiple input reference panels

## [ld_heatmap.py](ld_heatmap.py)
Creates a heatmap of the LD pattern among the variants in a genotype matrix.

## [manhattan.py](manhattan.py)
Create a manhattan plot to visualize the results of a GWAS.

## [plot_phenotypes.py](plot_phenotypes.py)
Performs a linear regression and outputs a plot for phenotypes from `generate_phenotypes.py`.

## [merge_plink.py](merge_plink.py)
Merge variants from two PGEN files that have the same set of samples.

## [run_FINEMAP.R](run_FINEMAP.R)
Run FINEMAP with a set of genotypes and phenotypes. Compute summary statistics, then call the `finemap_1p4.R` script.

## [run_SuSiE.R](run_SuSiE.R)
Run SuSiE with a set of genotypes and phenotypes. Follow the same interface as `run_FINEMAP.R`.

## [summarize_results.R](summarize_results.R)
This R script summarizes the output of several executions of FINEMAP and SuSiE by producing several useful visualizations.

## [utils.R](utils.R)
Define some useful methods for reading various file types into R. This script is sourced by other R scripts but never executed directly.
