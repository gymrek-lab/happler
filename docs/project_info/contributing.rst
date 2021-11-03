.. _project_info-contributing:

============
Contributing
============

Contributions are welcome and greatly appreciated!


------------
Dev Setup
------------

Follow these steps to set up a development environment.

1. Create a conda environment with ``poetry`` and ``sphinx``

    .. code-block:: console

        conda create -n happler-dev -c conda-forge poetry sphinx
2. Activate the environment

    .. code-block:: console

        conda activate happler-dev
3. Install our dependencies

    .. code-block:: console

        poetry install

---------------------
Managing Dependencies
---------------------
Run ``poetry help`` to read about the suite of commands it offers for managing dependencies.

For example, to add a pypi dependency to our list and install it, just run

    .. code-block:: console

        poetry add <dependency>

-----------
Code Checks
-----------
Before creating your pull request, please do the following.

1. Format the code correctly

    .. code-block:: console

        poetry run black .

2. If you made changes to the docs, check that they appear correctly.

    .. code-block:: console

        make -C docs html
        open docs/_build/html/index.html

3. Run all of the tests

    .. code-block:: console

        poetry run pytest tests/
