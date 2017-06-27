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

    This will clone both the FALCON-integrate and FALCON_unzip repositories and build a virtualenv including all FALCON daligner/dazzler dependencies.


.. IMPORTANT::


    FALCON_unzip also depends on several commands that are bundled with `SMRTLink <http://www.pacb.com/support/software-downloads/>`_
    and `MUMmer <https://sourceforge.net/projects/mummer/>`_. The easiest way to ensure these dependencies are in place
    and FALCON_unzip will run as expected is to

        1. Download and install `MUMmer 3.23 <https://sourceforge.net/projects/mummer/>`_
        2. Download and install `SMRTLink 4.0.0 <http://www.pacb.com/support/software-downloads/>`_


        .. code-block:: bash

            $ wget https://downloads.pacbcloud.com/public/software/installers/smrtlink_4.0.0.190159.zip
            $ unzip smrtlink_4.0.0.190159.zip   ### the password is SmrT3chN
            $ ./smrtlink_4.0.0.190159.run --rootdir smrtlink --smrttools-only

    If you don't care about unzipping / phasing your genome, or you are working with a haploid organism, this step is unnecessary.



Testing
-------

If you're reading this section, we are assuming that you have successfully completed the above steps and would now like
to test your install before you launch a FALCON/FALCON_unzip job on a large genome.

First, let's source our brand new virtualenv

.. code-block:: bash

    $ source /path/to/your/install/dir/fc_env/bin/activate

Second, let's ensure our small testing dataset is in place

.. code-block:: bash

    $ cd /path/to/your/install/dir/src/FALCON-integrate/FALCON-examples
    $ ../git-sym/git-sym update run/greg200k-sv2   ### ensure testdata in place

Third, we need to ensure the unarchived smrttools bin directory ``smrtlink/smrtcmds/bin`` (from instructions above) is specified correctly in your :ref:`fc_unzip.cfg`

.. code-block:: bash

    $ cd FALCON_examples/run/greg200k-sv2
    $ sed -i 's|^smrt_bin=.*!|smrt_bin=/path/to/smrtlink/smrtcmds/bin|g' fc_unzip.cfg

Fourth - we need to ensure the mummer & samtools (part of smrttools) binaries are available in your ``$PATH``

.. code-block:: bash

    $ export PATH=/path/to/mummer/bindir:$PATH
    $ export PATH=$PATH:/path/to/your/smrtlink/smrtcmds/bin

.. IMPORTANT::

    if you specify PATH=/path/to/smrtlink/smrtcmds/bin:$PATH instead of what's above you will run in to python import problems.

Fifth - let's test!

.. code-block:: bash

    $ fc_run fc_run.cfg
    $ fc_unzip.py fc_unzip.cfg
    $ fc_quiver.py fc_unzip.cfg

If you don't see any errors, you will have successfully assembled, unzipped, and polished a small test dataset. At this
point you should be ready to confidently launch a larger genome assembly.
