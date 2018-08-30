.. _hgap4_adapt:

.. caution:: These documents refer to an obsolete way of installing and running FALCON. They will remain up for historical context and for individuals still using the older version of FALCON/FALCON_unzip.

.. attention:: The current PacBio Assembly suite documentation which includes new bioconda instructions for installing FALCON, FALCON_unzip and their associated dependencies can be found here `pb_assembly <http://github.com/gconcepcion/pb-assembly>`_


============
Running Unzip on HGAP4 output
============

Overview
--------

``HGAP4`` is a ``FALCON``-based assembly pipeline, available through the SMRT Link interface. The pipeline itself encapsulates *de novo* assembly and polishing of the resulting contigs, but *not* the ``FALCON-unzip`` process as well. ``FALCON-unzip`` is currently available as a standalone tool, runnable only via command line.

Although ``HGAP4`` runs ``FALCON`` under the hood, the folder structure it generates is different than that of ``FALCON``. The ``FALCON-unzip``, however, requires the assembly folders to be formatted in the ``FALCON``-style.

This tutorial describes the necessary steps required to adjust the HGAP4 output to be compatible with a form required by ``FALCON-unzip``.


In brief, the majority of work required to adjust the ``HGAP4`` output to a ``FALCON``-compatible directory structure is implemented in a script called ``hgap4_adapt``. This script lives in the ``FALCON`` repository.

The complete process is composed of the following steps:

1. Installing ``FALCON`` and ``FALCON-unzip``.
2. Running ``hgap4_adapt``.
3. Creating the ``fc_unzip.cfg`` configuration file for ``FALCON-unzip``.
4. Creating the ``input.fofn`` and ``input_bam.fofn``.
5. Running ``FALCON-unzip``.

**IMPORTANT:** ``FALCON-unzip`` can only be run on HGAP4 jobs which had the ``Save Output for Unzip`` option turned on. It is not possible to run ``FALCON-unzip`` otherwise, because critical files will be missing from your job's output.


1. Installing ``FALCON`` and ``FALCON-unzip``
---------------------------------------------

The latest versions of ``FALCON`` and ``FALCON-unzip`` are available as precompiled Linux binaries. The easiest approach to installing them is through a wrapper script, described here:

:ref:`Quick Start<quick_start>`

Follow this approach to set-up the environment before moving on to step 2.

Alternatively, one can install the binaries manually by following the instructions here:
https://github.com/PacificBiosciences/FALCON_unzip/wiki/Binaries


2. Running ``hgap4_adapt``
--------------------------

Once the ``FALCON`` installation was successful, one needs to activate the installation environment to make the ``hgap4_adapt`` script available. This will also activate  ``FALCON`` and ``FALCON-unzip``. To verify the installation, run the following:

.. code-block:: bash

    source /path/to/your/install/dir/fc_env/bin/activate
    python -m falcon_kit.mains.hgap4_adapt --help

..

If everything was successful, this should output verbose usage information to screen. After this is set-up and working, adapting an existing HGAP4 run is as simple as the following example (take note of the dummy path, and replace it with a real one):

.. code-block:: bash

    source /path/to/your/install/dir/fc_env/bin/activate

    job_dir=/path/to/your/hgap4/job/123/123456/
    mkdir –p example1
    cd example1
    python -m falcon_kit.mains.hgap4_adapt --job-output-dir=${job_dir}

..

The result should be visible in the ``example1`` directory - it should now be populated to folders resembling a typical ``FALCON`` assembly run.


3. Creating the ``fc_unzip.cfg`` configuration file for ``FALCON-unzip``
------------------------------------------------------------------------

For help on .cfg files, please take a look at these Wiki pages:

- https://github.com/PacificBiosciences/FALCON/wiki
- https://github.com/PacificBiosciences/FALCON_unzip/wiki
- https://github.com/PacificBiosciences/FALCON-integrate/wiki


4. Creating the ``input.fofn`` and ``input_bam.fofn``
-----------------------------------------------------

The ``input.fofn`` file ("file of file names") contains the paths to files containing plain FASTA sequences of your raw reads, one file per row. All raw reads in the FASTA format should be available in your job dir:

.. code-block:: bash

    job_dir=/path/to/your/hgap4/job/123/123456/
    mkdir –p example1
    cd example1

    echo "${job_dir}/tasks/pbcoretools.tasks.gather_fasta-1/file.fasta" > input.fofn

..


The ``input_bam.fofn`` is required for the polishing step. This file is composed of a list of all BAM files from the input dataset which was provided to the initial HGAP4 run:

.. code-block:: bash

    source /path/to/your/install/dir/fc_env/bin/activate

    job_dir=/path/to/your/hgap4/job/123/123456/
    mkdir –p example1
    cd example1

    dataset summarize ${job_dir}/tasks/pbcoretools.tasks.filterdataset-0/filtered.subreadset.xml | grep -E "*.bam$" > input_bam.fofn

..

5. Running ``FALCON-unzip``
---------------------------

Before running ``FALCON-unzip``, the adapted folder structure should be similar to the following:

.. code-block:: bash

    $ cd example1
    $ ls | xargs -n 1
        0-rawreads
        1-preads_ovl
        2-asm-falcon
        fc_unzip.cfg
        input_bam.fofn
        input.fofn

..

Finally, to run ``FALCON-unzip``, do the following:

.. code-block:: bash

    source /path/to/your/install/dir/fc_env/bin/activate
    cd example1
    fc_unzip.py fc_unzip.cfg
    fc_quiver.py fc_unzip.cfg

..
