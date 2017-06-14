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

FALCON_unzip also depends on several commands that are bundled with SMRTAnalysis 4.0.0+.

They can be downloaded, built and installed manually (blasr/samtools/GenomicConsensus)

However, by far the easiest method is to download and install SMRTAnalysis tools. Once unarchived, you will use the fully
resolved PATH to the unarchived BIN directory ``/path/to/smrttools/install/smrtcmds/bin`` as the ``smrt_bin`` parameter in your FALCON_unzip
configuration file :ref:`fc_unzip.cfg`

.. code-block:: bash

    $ wget https://downloads.pacbcloud.com/public/software/installers/smrtlink_4.0.0.190159.zip
    $ unzip smrtlink_4.0.0.190159.zip   ### the password is SmrT3chN
    $ ./smrtlink_4.0.0.190159.run --rootdir smrtlink --smrttools-only

The ``--smrttools-only`` flag will allow you to install just the SMRTanalysis tools, skipping the full analysis suite
 configuration. Refer to this `PDF <http://programs.pacificbiosciences.com/e/1652/e-Installation--v4-0-0--v2-pdf/3rvmzg/507864561>`_
 if you would like to configure a fully functional SMRTLink 4.0.0 install
