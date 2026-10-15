.. _formats-haplotypes:


Haplotypes
==========

We utilize v0.2.0 of the `haptools .hap file format <https://haptools.readthedocs.io/en/stable/formats/haplotypes.html>`_ for haplotype data.

There are several additional "extra" attributes that are included in our haplotype files, which are described below. These attributes are not part of the official haptools format specification, but are included as helpful information.

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
     - Effect size
     - float
     - The effect size of the association between this haplotype and the trait
   * - 2
     - p-value
     - float
     - The p-value of the association between this haplotype and the trait
   * - 3
     - score
     - float
     - The delta BIC score of the final variant when it was added to this haplotype

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
     - Score
     - float
     - The delta BIC of adding this variant to the haplotype, up until this point
