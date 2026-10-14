.. _manual-main:

happler
=======

Happler is a haplotype-based fine-mapping method. This means that it tests for associations between a trait and haplotypes (ie sets of correlated SNPs instead of individual SNPs).

Commands
~~~~~~~~

* :doc:`happler run </commands/run>`: Run the haplotype-based fine-mapping method

* :doc:`happler transform </commands/transform>`: Transform a genotype matrix via a haplotype


Detailed information about each command can be found in the *Commands* section of our documentation. Examples there utilize a set of example files described :doc:`here </project_info/example_files>`.

Logging
~~~~~~~

All commands output log messages to standard error. The universal ``--verbosity`` flag controls the level of detail in our logging messages. By default, this is set to ``INFO``, which will yield errors, warnings, and info messages. To get more detailed messages, set it to ``DEBUG``. To get only error messages, set it to ``ERROR``. To get errors *and* warnings, set it to ``WARNING``. Refer to `the Python documentation on logging levels <https://docs.python.org/3/library/logging.html#levels>`_ for more information.

Contributing
~~~~~~~~~~~~

We gladly welcome any contributions to ``happler``!

Please read our :doc:`contribution guidelines </project_info/contributing>` and then submit a `Github issue <https://github.com/gymrek-lab/happler/issues>`_.

.. toctree::
   :caption: Inputs
   :name: inputs
   :hidden:
   :maxdepth: 1

   executing/inputs.rst

.. toctree::
   :caption: Outputs
   :name: outputs
   :hidden:
   :maxdepth: 1

   executing/outputs.rst

.. toctree::
   :caption: API
   :name: api
   :hidden:
   :maxdepth: 1

   api/modules

.. toctree::
   :caption: Project Info
   :name: project-info
   :hidden:
   :maxdepth: 1

   project_info/contributing
