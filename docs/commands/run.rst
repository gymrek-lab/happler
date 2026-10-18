.. _commands-run:


run
===

Run the haplotype construction step.

Overview
~~~~~~~~

The run command builds a forest of multiple trees from the input genotypes and
phenotypes. Each tree produces haplotypes, which are written to a ``.hap``
file for downstream analysis.

At each locus, happler uses phased genotypes to identify trait-associated
haplotypes: sets of alleles at one or more variants that occur in cis. Alleles
are selected greedily to strengthen the association between the haplotype and
the target trait.

Candidate haplotypes are identified with a binary tree of
linear regression models. Each node represents a variant, and its branches
represent the REF and ALT alleles. A path from the root to a node
defines the haplotype at that node; a haplotype with only one variant is
equivalent to an individual variant.

The root variant is chosen as the variant with the strongest trait association,
based on the Bayesian Information Criterion (BIC). To extend a haplotype,
happler evaluates candidate alleles at variants not already in the haplotype
and chooses the extension with the best association. Each candidate is encoded
as a haplotype matrix containing the number of copies (0, 1, or 2) carried by
each individual. Tree building stops when an extension does not improve model
fit enough according to the configured BIC thresholds, or when it no longer
provides sufficient evidence against independent effects of the parent
haplotype and candidate variant.

To account for multiple independent signals, happler builds multiple trees and
regresses the effects of haplotypes found by other trees out of the phenotype
before building each tree. The resulting output is the union of the haplotypes
from all trees in the forest. Final output can exclude haplotypes below the configured
association threshold (``--out-thresh``) and haplotypes consisting of a single variant (``--remove-SNPs``).

Usage
~~~~~

.. code-block:: bash

	happler run \
	--pheno TEXT \
	--pheno-name \
	--region TEXT \
	--sample SAMPLE --sample SAMPLE \
	--samples-file FILENAME \
	--discard-multiallelic \
	--discard-missing \
	--chunk-size INT \
	--maf FLOAT \
	--hap-maf FLOAT \
	--max-signals INT \
	--max-iterations INT \
	--ld-prune-thresh FLOAT \
	--out-thresh FLOAT \
	--show-tree \
	--remove-SNPs \
	--output PATH \
	--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
	GENOTYPES PHENOTYPES

Input
~~~~~

Phased genotypes may be in VCF/BCF or PLINK2 PGEN formats. Phenotypes
must be a tab-separated phenotype file containing sample IDs and phenotype
values. Refer to the :doc:`genotype format documentation </formats/genotypes>` and :doc:`phenotype format documentation </formats/phenotypes>` for more information.

By default, the first phenotype column is used. Use ``--pheno`` to
select a phenotype by index, or combine it with ``--pheno-name`` to select a
column by name.

The input can be restricted to a genomic region with ``--region`` or to
selected samples with ``--sample`` or ``--samples-file``. ``--sample`` may be
repeated, but it cannot be combined with ``--samples-file``.

Input selection
~~~~~~~~~~~~~~~

``--pheno TEXT``
	Select the phenotype by its zero-based column index. The default is
	``0``, the first phenotype column.

``--pheno-name``
	Treat the value passed to ``--pheno`` as a phenotype column name instead
	of an index.

``--region TEXT``
	Restrict genotypes to a genomic region, such as ``chr1:1234-34566`` or
	``chr7``. The VCF must be indexed and its sequence names must match.

``--sample SAMPLE``
	Restrict the analysis to the named sample. This option may be repeated.

``--samples-file FILENAME``
	Read sample IDs from a file, one per line. This cannot be combined with
	``--sample``.

Output
~~~~~~

The command outputs haplotypes in the ``.hap`` format to standard output by default.

Use ``--output`` to write to a file. Use ``--show-tree`` to additionally
write the constructed forest in Graphviz DOT format. When ``--output`` names a
file, the DOT output uses the same path with a ``.dot`` suffix.

Filtering and memory
~~~~~~~~~~~~~~~~~~~~

``--discard-multiallelic``
	Discard multiallelic variants instead of reporting them as an error.

``--discard-missing``
	Discard samples with missing genotypes for required variants.

``--maf FLOAT``
	Discard input variants whose minor allele frequency is below ``FLOAT``.

``--hap-maf FLOAT``
	Only use variants with a minor allele frequency above ``FLOAT`` when
	building haplotypes. By default, use the MAF equivalent to a minor allele
	count of 20.

``--chunk-size INT``
	Read PGEN variants in chunks of ``INT`` variants. Chunking reduces memory
	use at the cost of additional runtime.

Forest and output options
~~~~~~~~~~~~~~~~~~~~~~~~~

``--max-signals INT``
	Build a forest with up to ``INT`` trees, corresponding to the maximum
	number of expected causal signals.

``--max-iterations INT``
	Repeat tree building ``INT`` times.

``--ld-prune-thresh FLOAT``
	Prune leaf nodes based on linkage disequilibrium with their siblings when
	the LD exceeds ``FLOAT``.

``--out-thresh FLOAT``
	Only output haplotypes meeting this single-SNP significance threshold.

``--show-tree``
	Write the forest as a Graphviz DOT file in addition to the haplotype output.

``--remove-SNPs``
	Do not output haplotypes containing only a single variant. This also omits
	single-variant trees from DOT output.

``--output PATH``
	Write haplotypes to ``PATH`` instead of standard output.

``--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET]``
	Set the logging verbosity. The default is ``INFO``.

Examples
~~~~~~~~

Run haplotype construction and write haplotypes to standard output:

.. code-block:: bash

	happler run \
		tests/data/simple.vcf \
		tests/data/simple.tsv \
		> simple.hap

Write haplotypes to a file and use three signals and three iterations:

.. code-block:: bash

	happler run \
		--max-signals 3 \
		--max-iterations 3 \
		--output haplotypes.hap \
		tests/data/simple.vcf \
		tests/data/simple.tsv

Select a phenotype by column name and restrict the analysis to a region:

.. code-block:: bash

	happler run \
		--pheno-name \
		--pheno hap \
		--region 19:45401409-46401409 \
		--output haplotypes.hap \
		tests/data/19_45401409-46401409_1000G.pgen \
		tests/data/19_45401409-46401409_1000G.pheno

Write the forest as a DOT file and omit single-variant haplotypes:

.. code-block:: bash

	happler run \
		--show-tree \
		--remove-SNPs \
		--output haplotypes.hap \
		tests/data/simple.vcf \
		tests/data/simple.tsv

All files used in these examples are described :doc:`here </project_info/example_files>`.

Detailed Usage
~~~~~~~~~~~~~~

.. click:: happler.__main__:main
   :prog: happler
   :nested: full
   :commands: run
