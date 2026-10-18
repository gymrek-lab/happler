.. _commands-transform:


transform
=========

Transform a genotype matrix using a haplotype from a ``.hap`` file.

The command removes the variants that define the selected haplotype and
returns pseudo-genotypes for that haplotype extended by each retained variant.
Each output variant therefore represents the number of chromosomes carrying
the selected haplotype and the requested allele at that variant. This is
useful for testing haplotype extensions in downstream association or
conditional analyses.

One haplotype is transformed at a time. If ``--hap-id`` is not specified, the
first haplotype in the input file is used.

Usage
~~~~~

.. code-block:: bash

	happler transform \
	--region TEXT \
	--sample SAMPLE --sample SAMPLE \
	--samples-file FILENAME \
	--var-id VAR_ID --var-id VAR_ID \
	--discard-multiallelic \
	--discard-missing \
	--chunk-size INT \
	--maf FLOAT \
	--hap-id HAP_ID \
	--allele INT \
	--output PATH \
	--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
	GENOTYPES HAPLOTYPES

Input
~~~~~

Phased genotypes may be in VCF/BCF or PLINK2 PGEN formats. Haplotypes
must be a file in the :doc:`.hap format </formats/haplotypes>`.

The genotype and haplotype files must contain compatible variant IDs. The
variants defining the selected haplotype are included automatically. Use
``--var-id`` to retain additional variants in the transformed output.

The input can be restricted to a genomic region with ``--region`` or to
selected samples with ``--sample`` or ``--samples-file``. ``--sample`` may be
repeated, but it cannot be combined with ``--samples-file``.

Output
~~~~~~

The command outputs pseudo-genotypes to standard output in VCF format by
default. Use ``--output`` to write to a file. The output format is selected
from the filename suffix (``.vcf``, ``.bcf`` or ``.pgen``). See :ref:`the genotype format
documentation <formats-genotypes>` for more information.

Filtering and memory
~~~~~~~~~~~~~~~~~~~~

``--discard-multiallelic``
	Discard multiallelic variants instead of reporting them as an error.

``--discard-missing``
	Discard samples with missing genotypes for required variants.

``--maf FLOAT``
	Discard variants whose minor allele frequency is below ``FLOAT``. The
	filter is applied to the input and transformed output.

``--chunk-size INT``
	Read PGEN variants in chunks of ``INT`` variants. Chunking reduces memory
	use at the cost of additional runtime.

Transformation options
~~~~~~~~~~~~~~~~~~~~~~

``--hap-id HAP_ID``
	Select a haplotype by ID. If omitted, the first haplotype is selected.

``--allele INT``
	Select the allele of the additional causal variant to append to the
	haplotype. The default is ``0``; use ``1`` for the alternate allele.

Examples
~~~~~~~~

Transform the first haplotype in a PGEN file and write VCF to standard output:

.. code-block:: bash

	happler transform \
		tests/data/19_45401409-46401409_1000G.pgen \
		tests/data/19_45401409-46401409_1000G.hap \
		> transformed.vcf

Select a particular haplotype and write the result to a VCF file:

.. code-block:: bash

	happler transform \
		--hap-id H0 \
		--output transformed.vcf \
		tests/data/19_45401409-46401409_1000G.pgen \
		tests/data/19_45401409-46401409_1000G.hap

Append the alternate allele (instead of the default REF allele) while transforming:

.. code-block:: bash

	happler transform \
		--hap-id H0 \
		--allele 1 \
		--output transformed.vcf \
		tests/data/19_45401409-46401409_1000G.pgen \
		tests/data/19_45401409-46401409_1000G.hap

Restrict the transformation to a region and selected samples:

.. code-block:: bash

	happler transform \
		--region 19:45401409-46401409 \
		--sample HG00096 \
		--sample HG00097 \
		--output transformed.vcf \
		tests/data/19_45401409-46401409_1000G.pgen \
		tests/data/19_45401409-46401409_1000G.hap

Read samples from a file and produce PGEN output:

.. code-block:: bash

	happler transform \
		--samples-file samples.txt \
		--chunk-size 1000 \
		--output transformed.pgen \
		tests/data/19_45401409-46401409_1000G.pgen \
		tests/data/19_45401409-46401409_1000G.hap

Retain an additional variant by ID:

.. code-block:: bash

	happler transform \
		--hap-id H0 \
		--var-id rs8106922 \
		--output transformed.vcf \
		tests/data/19_45401409-46401409_1000G.pgen \
		tests/data/19_45401409-46401409_1000G.hap

All files used in these examples are described :doc:`here </project_info/example_files>`.

Detailed Usage
~~~~~~~~~~~~~~

.. click:: happler.__main__:main
   :prog: happler
   :nested: full
   :commands: transform
