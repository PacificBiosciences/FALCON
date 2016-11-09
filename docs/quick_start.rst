.. _quick_start:

Quick Start Guide
=================


FALCON
------

Installation
^^^^^^^^^^^^

By far, the quickest way to get up and running with FALCON is to ensure you have python 2.7.x installed and in your
path, then follow these simple instructions to clone and install the FALCON-integrate_ github repository:

.. code-block:: bash

    $ git clone git://github.com/PacificBiosciences/FALCON-integrate.git
    $ cd FALCON-integrate
    $ make init
    $ source env.sh
    $ make config-edit-user
    $ make -j all
    $ make test  # to run a simple one


.. _FALCON-integrate: https://github.com/PacificBiosciences/FALCON-integrate

Usage
^^^^^

To run a job with your newly installed FALCON-integrate package, simply add the ``fc_env/bin`` directory to your
``$PATH`` and call invoke :doc:`fc_run.py` with your ``fc_run.cfg`` as the argument.

.. code-block:: bash

    falcon_jobdir$ export PATH=/path/to/FALCON-integrate/fc_env/bin:$PATH
    falcon_jobdir$ fc_run.py fc_run.cfg


FALCON_unzip
------------

Installation
^^^^^^^^^^^^

Installation is very simple, as you are going to install FALCON_unzip into your previously created FALCON ``fc_env``.

.. code-block:: bash

    $ git clone https://github.com/PacificBiosciences/FALCON_unzip.git
    $ export PATH=/path/to/FALCON-integrate/fc_env/bin:$PATH
    $ cd FALCON_unzip
    $ pip install --user ./


Usage
~~~~~

FALCON_unzip operates on a completed FALCON job directory. All that's needed is to create an :ref:`fc_unzip.cfg` and place
it at the ``$ROOT`` of your FALCON job. Then invoke with :ref:`fc_unzip.py`:

.. code-block:: bash

    falcon_jobdir$ export PATH=/path/to/FALCON-integrate/fc_env/bin:$PATH
    falcon_jobdir$ fc_unzip.py fc_unzip.cfg

*Note: FALCON_unzip also depends on several commands that are bundled with SMRTAnalysis 3.0.x, not currently
available for the public. The workaround is to download each tool independently and link the binaries to a path that
you set as in :ref:`fc_unzip.cfg` as your ``smrt_bin`` parameter.

1. `samtools <https://github.com/samtools/samtools>`_
2. `pbalign <https://github.com/PacificBiosciences/pbalign>`_
3. `pbbam <https://github.com/PacificBiosciences/pbbam>`_
4. `variantCaller <https://github.com/PacificBiosciences/GenomicConsensus>`_

