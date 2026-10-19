.. _project_info-installation:

============
Installation
============

Using pip
---------

You can install ``happler-fm`` from PyPI using ``pip``. This provides the
``happler`` Python package and ``happler`` command.

.. note::
   We recommend using ``pip >= 20.3``.

   .. code-block:: bash

      pip install 'pip>=20.3'

.. code-block:: bash

   pip install happler-fm

Using conda
-----------

We also support installing ``happler`` from bioconda using ``conda``.

.. code-block:: bash

   conda install -c conda-forge -c bioconda happler

Installing the latest, unreleased version
-----------------------------------------
Can't wait for us to tag and release our most recent updates? You can install ``happler`` directly from the ``main`` branch of our Github repository using ``pip``.

.. code-block:: bash

   pip install --upgrade --force-reinstall git+https://github.com/gymrek-lab/happler.git

Development installation
------------------------
If you would like to edit and work with the source code, please refer to :ref:`these instructions <dev-setup-instructions>` for a different installation setup.
