.. _quick_start:

Quick Start Guide
=================


FALCON
------

Installation
^^^^^^^^^^^^

The quickest way to install FALCON + FALCON_unzip is to download and run this install script:

:down`load:`install_unzip.sh <scripts/install_unzip.sh>`

.. code-block:: bash

    $ bash -ex install_unzip.sh /path/to/your/install/dir

This will clone both the FALCON-integrate and FALCON_unzip repositories and build a virtualenv including all FALCON daligner/dazzler dependencies.

Install script dependencies (Ubuntu/CentOS):

 1. python/2.7.x
 2. virtualenv 13.0.1 (but should work with older)
 3. git binaries in your $PATH


Additional FALCON_unzip dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FALCON_unzip also depends on several commands that are bundled with SMRTAnalysis 3.0.x, but not currently
available for the public. The workaround is to download each tool independently and link the binaries to a path that
you set in :ref:`fc_unzip.cfg` as your ``smrt_bin`` parameter. This is not necessary if you are not interested in
phasing your genome.

1. `samtools <https://github.com/samtools/samtools>`_
2. `pbalign <https://github.com/PacificBiosciences/pbalign>`_
3. `GenomicConsensus v1.1.0@654d0276d4a03f269cd1a14ddd7dfd0f54bede45 <https://github.com/PacificBiosciences/GenomicConsensus/tree/654d0276d4a03f269cd1a14ddd7dfd0f54bede45>`_

** Note: the script ``makePbi.py`` is only available from GenomicConsensus 1.1.0 and earlier. To sync the branch
that has this script, follow these steps:

.. code-block:: bash

    $ git clone https://github.com/PacificBiosciences/GenomicConsensus.git
    $ cd GenomicConsensus
    $ git reset --hard 654d0276d4a03f269cd1a14ddd7dfd0f54bede45
