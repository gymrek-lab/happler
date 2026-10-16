.. _commands-transform:


transform
=========

Append a new allele to all haplotypes in a genotypes file. Return the new set of genotypes.

Usage
~~~~~
.. code-block:: bash

	happler transform \
	--region TEXT \
	--sample SAMPLE --sample SAMPLE \
	--samples-file FILENAME \
	--chunk-size INT \
	--discard-missing \
	--output PATH \
	--verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
	GENOTYPES HAPLOTYPES

Input
~~~~~
Genotypes must be specified in VCF and haplotypes must be specified in the :doc:`.hap file format </formats/haplotypes>`.

Alternatively, you may specify genotypes in PLINK2 PGEN format. Just use the appropriate ".pgen" file extension in the input. See the documentation for genotypes in :ref:`the format docs <formats-genotypesplink>` for more information.

Output
~~~~~~
Transform outputs *psuedo-genotypes* in VCF, but you may request genotypes in PLINK2 PGEN format, instead. Just use the appropriate ".pgen" file extension in the output path. See the documentation for genotypes in :ref:`the format docs <formats-genotypesplink>` for more information.

Examples
~~~~~~~~
# TODO: Add more examples of usage, including with PGEN files.

All files used in these examples are described :doc:`here </project_info/example_files>`.


Detailed Usage
~~~~~~~~~~~~~~

.. click:: happler.__main__:main
   :prog: happler
   :nested: full
   :commands: transform