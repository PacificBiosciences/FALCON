.. _quick_start:

Quick Start Guide
=================


FALCON
------

Installation
^^^^^^^^^^^^

The quickest way to install FALCON + FALCON_unzip is to download and run this install script:

:download:`install_unzip.sh <scripts/install_unzip.sh>`

.. code-block:: bash

    $ bash -ex install_unzip.sh /path/to/your/install/dir

This will clone both the FALCON-integrate and FALCON_unzip repositories and build a virtualenv including all FALCON daligner/dazzler dependencies.

Install script dependencies (Ubuntu/CentOS):

 1. python/2.7.x
 2. virtualenv 13.0.1 (but should work with older)
 3. git binaries in your $PATH


Additional FALCON_unzip dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FALCON_unzip also depends on several commands that are bundled with SMRTAnalysis 4.0.0+. You can download/build these
dependencies manually, (blasr/samtools/pbalign/GenomicConsensus) or by far the easiest option is to download the
SMRTAnalysis 4.0.0 zip file, unarchive, deploy and use the full PATH to the `smrtcmds/bin` directory as your ``smrt_bin``
parameter in :ref:`fc_unzip.cfg`. If you are not interested in unzipping and phasing your genome, installing these
additional FALCON_unzip dependencies is unnecessary.

.. code-block:: bash

    $ wget https://downloads.pacbcloud.com/public/software/installers/smrtlink_4.0.0.190159.zip
    $ unzip smrtlink_4.0.0_190159.zip   ### the password is SmrT3chN
    $ ./smrtlink_4.0.0.190159.run

For a simple local install being used exclusively as binaries for FALCON, you can select default options across the board.

If you want to install a fully functional SMRTAnalysis GUI suite, you should refer to the Release Documentation section
`here <http://www.pacb.com/support/software-downloads/>`_.

TLDR: pdf of the SMRTAnalysis installation process `here <http://programs.pacificbiosciences.com/e/1652/e-Installation--v4-0-0--v2-pdf/3rvmzg/507864561>`_
