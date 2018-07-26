.. _quick_start:

Quick Start Guide
=================

Installation
------------

The quickest way to install FALCON + FALCON_unzip is to download and run this install script:


============================================================== ====================== =======
Install script                                                      Tarball date      Status
============================================================== ====================== =======
:download:`install_unzip.sh <scripts/install_unzip_180312.sh>` :ref:`3122018tarball`
:download:`install_unzip.sh <scripts/install_unzip_180725.sh>` :ref:`7252018tarball`   Beta
============================================================== ====================== =======


.. code-block:: bash

    $ bash -ex install_unzip.sh /path/to/your/install/dir


.. NOTE::

    This will clone the FALCON-integrate repository, FALCON_unzip binaries, build a virtualenv and launch a small test case assembly to ensure successful installation.

    samtools_ must be installed separately and in your $PATH.

    If you don't see any errors, you will have installed FALCON/FALCON_unzip and successfully assembled and unzipped a small test dataset. At this point you should be ready to confidently launch a larger genome assembly.

To activate your FALCON_unzip virtualenv in the future:

.. code-block:: bash

    $ source /path/to/your/install/dir/fc_env/bin/activate


.. _samtools: http://www.htslib.org/download/
.. _minimap2: https://github.com/lh3/minimap2
