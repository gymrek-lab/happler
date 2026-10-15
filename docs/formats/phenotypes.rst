.. _formats-phenotypes:


Phenotypes and Covariates
=========================

Phenotype file format
---------------------
Phenotypes are expected to follow `the PLINK2 .pheno file format <https://www.cog-genomics.org/plink/2.0/input#pheno>`_. This is a
tab-separated format where the first column corresponds to the sample ID, and subsequent columns contain each of your phenotypes.

The first line of the file corresponds with the header and must begin with ``#IID``.
The names of each of your phenotypes belong in the subsequent columns of the header.

See `tests/data/19_45401409-46401409_1000G.pheno <https://github.com/gymrek-lab/happler/blob/main/tests/data/19_45401409-46401409_1000G.pheno>`_ for an example of a phenotype file:

.. include:: ../../tests/data/19_45401409-46401409_1000G.pheno
    :literal:
    :start-line: 1
    :end-line: 6
