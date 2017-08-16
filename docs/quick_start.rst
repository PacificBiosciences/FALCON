.. _quick_start:

Quick Start Guide
=================

Installation
------------

The quickest way to install FALCON + FALCON_unzip is to download and run this install script:

:download:`install_unzip.sh <scripts/install_unzip.sh>`

.. code-block:: bash

    $ bash -ex install_unzip.sh /path/to/your/install/dir


.. NOTE::

    This will clone the FALCON-integrate repository, FALCON_unzip binaries, build a virtualenv that includes all necessary dependencies.
and launch a small test case assembly to ensure successful installation.

If you don't see any errors, you will have installed FALCON/FALCON_unzip and successfully assembled and unzipped a small test dataset.
At this point you should be ready to confidently launch a larger genome assembly.

To activate your FALCON_unzip virtualenv in the future:

.. code-block:: bash

    $ source /path/to/your/install/dir/fc_env/bin/activate
