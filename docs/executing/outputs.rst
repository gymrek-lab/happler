.. _executing-outputs:


Outputs
=========

Haplotype file format
---------------------
Discovered haplotypes will be written in a special ``.haps`` file format.

This is a tab-separated file composed of different types of lines. The first field of each line is a single, uppercase character denoting the type of line.

.. list-table::
   :widths: 25 25
   :header-rows: 1

   * - Type
     - Description
   * - #
     - Comment
   * - M
     - Metadata
   * - H
     - Haplotype
   * - V
     - Variant

``#`` Comment line
~~~~~~~~~~~~~~~~~~
Comment lines begin with ``#`` and are ignored.

``M`` Metadata
~~~~~~~~~~~~~~
Describes attributes of the file, itself, including the version of the format spec.

``H`` Haplotype
~~~~~~~~~~~~~~~
Haplotypes contain the following attributes:

.. list-table::
   :widths: 25 25 25 50
   :header-rows: 1

   * - Column
     - Field
     - Type
     - Description
   * - 1
     - Haplotype ID
     - int
     - Identifies a haplotype within a tree; unique across other haplotypes in the tree
   * - 2
     - Tree ID
     - int
     - A unique identifier for the tree that this haplotype belongs to
   * - 3
     - Effect size
     - float
     - The effect size of the association between this haplotype and the trait
   * - 4
     - p-value
     - float
     - The p-value of the association between this haplotype and the trait
   * - 5
     - PIP
     - float
     - The posterior inclusion probability for the causality of this haplotype

``V`` Variant
~~~~~~~~~~~~~
Variant lines must follow the haplotype to which they belong. They contain the following attributes:

.. list-table::
   :widths: 25 25 25 50
   :header-rows: 1

   * - Column
     - Field
     - Type
     - Description
   * - 1
     - Variant ID
     - string
     - The unique ID for this variant, as defined in the genotypes file
   * - 2
     - Allele
     - bool
     - The allele of this variant within the haplotype
   * - 3
     - Score
     - float
     - The importance of including this variant within the haplotype
