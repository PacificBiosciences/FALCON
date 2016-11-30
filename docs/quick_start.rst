.. _quick_start:

Quick Start Guide
=================


FALCON
------

Installation
^^^^^^^^^^^^

By far, the quickest way to get up and running with FALCON is to ensure you have python 2.7.x installed and in your
path, then follow these simple instructions to clone and install the FALCON-integrate_ github repository. This includes
proper versions of all necessary dependencies. Alternatively, you can install the FALCON github package manually,
but you will need to do some legwork to make sure you have the right versions of all dependencies installed.

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
``$PATH`` and invoke :doc:`fc_run.py` with your ``fc_run.cfg`` as the argument.

.. code-block:: bash

    falcon_jobdir$ export PYTHONUSERBASE=/path/to/FALCON-integrate/fc_env
    falcon_jobdir$ export PATH=$PYTHONUSERBASE/bin:$PATH
    falcon_jobdir$ fc_run.py fc_run.cfg


FALCON_unzip
------------

Installation
^^^^^^^^^^^^

Installation is very simple, as you are going to install FALCON_unzip into your previously created FALCON ``fc_env``.

.. code-block:: bash

    $ git clone https://github.com/PacificBiosciences/FALCON_unzip.git
    $ export PYTHONUSERBASE=/path/to/FALCON-integrate/fc_env
    $ cd FALCON_unzip
    $ python setup.py install --prefix=/path/to/FALCON-integrate/fc_env


Usage
^^^^^

FALCON_unzip operates on a completed FALCON job directory. All that's needed is to create an :ref:`fc_unzip.cfg` and place
it at the ``$ROOT`` of your FALCON job. Then invoke with :ref:`fc_unzip.py`:

.. code-block:: bash

    falcon_jobdir$ export PYTHONUSERBASE=/path/to/FALCON-integrate/fc_env
    falcon_jobdir$ export PATH=$PYTHONUSERBASE/bin:$PATH
    falcon_jobdir$ fc_unzip.py fc_unzip.cfg

`* Note: FALCON_unzip also depends on several commands that are bundled with SMRTAnalysis 3.0.x, not currently
available for the public. The workaround is to download each tool independently and link the binaries to a path that
you set in :ref:`fc_unzip.cfg` as your ``smrt_bin`` parameter.

1. `samtools <https://github.com/samtools/samtools>`_
2. `pbalign <https://github.com/PacificBiosciences/pbalign>`_
3. `GenomicConsensus v1.1.0@654d0276d4a03f269cd1a14ddd7dfd0f54bede45 <https://github.com/PacificBiosciences/GenomicConsensus/tree/654d0276d4a03f269cd1a14ddd7dfd0f54bede45>`_

`** Note: the script ``makePbi.py`` is only available from GenomicConsensus 1.1.0 and earlier. To sync the branch
that has this script, follow these steps:

.. code-block:: bash

    $ git clone https://github.com/PacificBiosciences/GenomicConsensus.git
    $ cd GenomicConsensus
    $ git reset --hard 654d0276d4a03f269cd1a14ddd7dfd0f54bede45
