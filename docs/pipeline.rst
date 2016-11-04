.. _pipeline:

.. image:: falcon_icon2.png
   :height: 200px
   :width: 200 px
   :alt: Falcon Assembler
   :align: right


FALCON Pipeline
===============


A FALCON job can be broken down into 3 steps:

1. Overlap detection and :term:`error correction` of rawreads
2. Overlap detection between corrected reads (:term:`pread`)
3. :term:`String Graph` assembly of corrected reads

Each step is performed in it's own subdirectory within the FALCON job

.. code-block:: bash

    falcon_job/
        ├── 0-rawreads     # Raw read error correction directory
        ├── 1-preads_ovl   # Corrected read overlap detection
        ├── 2-asm-falcon   # String Graph Assembly
        ├── scripts
        └── sge_log        # Job scheduler logs

The assembly process is driven by the script ``fc_run.py`` which should be run on a head node as it needs to persist
throughout the entire assembly process.
It takes as input a single :ref:`config <Configuration>` file typically named ``fc_run.cfg``, which references a list 
of fasta input files. The config file can be configured to
run locally, or submit to a job scheduler. However, if your dataset is anything larger than a bacterial sized
genome and unless you've tuned your system specifically for the organism you're trying to assemble, then most likely you
should be running on a cluster in order to more effectively leverage your computational resources.

The configuration file also allows you to control other aspects of your job such as how your compute resources are
distributed as well as set many parameters to help you reach an "optimized" assembly according to the nature of
your input data. Unfortunately at this point there is no "magic" way to auto-tune the parameters so you should
probably spend some time in the :ref:`configuration` section to understand what options are available to you.
Some example configuration files can be found :ref:`here <configuration>`



Step 1: Overlap detection and error correction of raw reads
-----------------------------------------------------------

The first step of the pipeline is to identify all overlaps in the raw reads. Currently this is all performed with
a modified version of Gene Myers' Daligner_.

In order to identify overlaps, your :term:`raw reads` must first be converted from fasta format into a Dazzler_
database. This is a very I/O intensive process and will be run from the node where ``fc_run.py`` was executed. If this
is an issue, you should submit the command with a wrapper script to your grid directly.

Once the database has been created and partitioned according to the parameters set in your
:ref:`fc_run.cfg <configuration>`, an all vs all comparison of the reads must be performed. Accordingly, due to the
all vs all nature of the search this is the most time consuming step in the assembly process. To walk through the
actual steps of this part of the pipeline you should take a look at `0-rawreads/prepare_rdb.sub.sh`.
Essentially it consists of running:

1. ``fasta2DB`` to format the database
2. ``DBSplit`` to partition the database
3. ``HPC.daligner`` to generate the ``daligner`` commands necessary for all-vs-all comparison

After overlaps have been detected, you will be left with many directories full of alignment files ``*.las`` containing
the information about the overlaps. After merging the alignment files, the next step is to error correct the reads
leveraging the overlap information. In the ``0-rawreads/preads`` directory you will find a series of scripts for
performing the error correction. The process basically consists of using ``LA4Falcon`` with a length cutoff and piping the
output to :ref:`fc_consensus.py <fc_consensus>` to generate a fasta file with corrected reads.

One can use :ref:`cns_concurrent_jobs` to control the maximum number of concurrent consensus jobs submitted to the job
management system. The ``out.XXXXX.fasta`` files produced are used as input for the next step in the pipeline.

.. code-block:: bash

    0-rawreads/
        ├── pre_assembly_stats.json   # pre-assembly stats
        ├── cns_done                  # consensus sentinel file
        ├── preads/                   # sub-dir for error correction
        ├── m_*/                      # dirs for all of the LA4Merge jobs
        ├── da_done                   # daligner sentinel file
        ├── job_*                     # dirs for all of the daligner jobs
        ├── length_cutoff             # text file with just the length cutoff
        ├── raw_reads.db              # dazzler DB file
        ├── rdb_build_done            # database construction sentinel file
        ├── run_jobs.sh               # listing of all overlap step commands
        ├── input.fofn                # list if your input *.fasta files
        ├── prepare_rdb.sh            # env wrapper script
        └── prepare_rdb.sub.sh        # driver script for this step in the pipeline



.. _Daligner: http://dazzlerblog.wordpress.com
.. _Dazzler: https://dazzlerblog.wordpress.com/2014/06/01/the-dazzler-db/


Step 2: Overlap detection of corrected reads
--------------------------------------------

Starting from the error corrected reads generated in the first step of the pipeline, we now need to perform an
additional overlap detection step. Depending on how well the error correction step proceeded as well as the how much
initial coverage was fed into the pipeline, the input data for this step should be significantly reduced at this
point. Thus, while still time consuming, the corrected read overlap detection step should proceed significantly faster.

The commands in this step of the pipeline are very similar to before albeit with different parameter settings to account
for the corrected nature of the :term:`pread`s. See ``1-preads_ovl/prepare_pdb.sub.sh`` for details on the parameters.

The only conceptual difference between the first and second overlap detection steps is that consensus calling is
only performed in the case of the initial raw read correction. After :term:`pread` overlap detection, it's simply a matter of
extracting the information from the corrected reads database ``DB2Falcon -U preads``.

.. code-block:: bash

    1-preads_ovl/
        ├── db2falcon_done          # sentinel file for final preads4falcon.fasta output
        ├── preads4falcon.fasta     # final corrected reads used in Assembly Graph
        ├── run_db2falcon.sh        # env wrapper script
        ├── run_db2falcon.sub.sh    # script to output preads from dazzler DB
        ├── p_merge_done            # sentinel file for *.las merging completion
        ├── m_*/                    # directories for LA4Merge jobs
        ├── da_done                 # sentinel file for completion of daligner jobs
        ├── job_*/                  # directories for daligner jobs
        ├── pdb_build_done          # sentinel file for pread DB building
        ├── preads.db               # preads dazzler DB
        ├── run_jobs.sh             # listing of all pread overlap job commands
        ├── prepare_pdb.sh          # env wrapper script
        ├── prepare_pdb.sub.sh      # driver script for this step of the pipeline
        └── input_preads.fofn       # list of your out.XXXXX.fasta's from previous step


Step 3: String Graph assembly
-----------------------------

The final step of the FALCON Assembly pipeline is generation of the final :term:`String Graph` assembly and output in
fasta format. There are 4 commands being run in the final phase of the FALCON assembly pipeline:

1. ``fc_ovlp_filter`` Filters overlaps based on the criteria you provided in fc_run.cfg
2. ``fc_ovlp_to_graph`` constructs an overlap graph of reads larger than the ``--min_len`` threshold provided
3. ``fc_graph_to_contig`` generates fasta files for contigs from the overlap graph.
4. ``fc_dedup_a_tigs`` removes duplicated associated contigs

You can see the details on the parameters used by inspecting ``2-asm_falcon/run_falcon_asm.sub.sh``
This step of the pipeline is very fast relative to the overlap detection steps. Sometimes it may be useful to run
several iterations of this step with different parameter settings in order to identify a "best" assembly.

The final output of this step is a fasta file of all of the primary contigs, ``p_ctg.fa`` as well as an associated contig
fasta file, ``a_ctg.fa`` that consists of all of the structural variants from the primary contig assembly.

.. code-block:: bash

    2-asm-falcon/
        ├── a_ctg_all.fa                 # all associated contigs, including duplicates
        ├── a_ctg_base.fa                #
        ├── a_ctg_base_tiling_path       #
        ├── a_ctg.fa                     # De-duplicated associated fasta file
        ├── a_ctg_tiling_path            # tiling path informaiton for each associated contig
        ├── falcon_asm_done              # FALCON Assembly sentinal file
        ├── p_ctg.fa                     # Fasta file of all primary contigs
        ├── p_ctg_tiling_path            # Tiling path of preads through each primary contig
        ├── c_path                       #
        ├── ctg_paths                    # corrected read paths for each contig
        ├── fc_ovlp_to_graph.log         # logfile for process of converting overlaps to assembly graph
        ├── utg_data                     #
        ├── sg_edges_list                # list of all edges
        ├── chimers_nodes                #
        ├── preads.ovl                   # List of all overlaps between preads
        ├── las.fofn                     # List of *.las files for input
        ├── run_falcon_asm.sh            # env wrapper script
        └── run_falcon_asm.sub.sh        # Assembly driver script



Supplementary Information
-------------------------
Supplemental command reference


Dazzler commands
----------------
These commands are part of Gene Meyer's Dazzler Suite of tools

.. _daligner:

daligner
++++++++
info

.. _DB2Falcon:

DB2Falcon
+++++++++
Used to dump dazzler preads.db into FASTA format for subsequent :term:`String Graph` assembly

.. _DB2Fasta:

DB2Fasta
++++++++
info

.. _DBdump:

DBdump
++++++
info

.. _DBdust:

DBdust
++++++

.. _DBsplit:

DBsplit
+++++++
The total number of jobs that are run is determined by how one "splits" the sequence database. You should read
Gene Myers's blog `Dazzler Blog <http://dazzlerblog.wordpress.com>` carefully to understand how the tuning options,
:ref:`pa_DBsplit_option` and :ref:`pa_HPCdaligner_option` work. Generally, for large genomes, you should use
``-s400`` (400Mb sequence per block) in :ref:`pa_DBsplit_option`. This will make a smaller number of jobs but each
job will run longer. However, if you have a job scheduler which limits how long a job can run, it might be
desirable to have a smaller number for the ``-s`` option.

.. _DBstats:
DBstats
+++++++

.. _fasta2DB:

fasta2DB
++++++++
info

.. _HPC.daligner:

HPC.daligner
++++++++++++
info

.. _LA4Falcon:

LA4Falcon
+++++++++
Output data from a Dazzler DB into fasta format for FALCON. You can supply the argument ``-H`` with an integer value
to filter reads below a given threshold.

.. _LAcheck:

LAcheck
+++++++

Check integrity of alignment files.

.. _LAmerge:

LAmerge
+++++++

The total number of jobs that are run is determined by how one "splits" the sequence database. You should read
Gene Myers's blog ( http://dazzlerblog.wordpress.com ) carefully to know how to tune the option pa_DBsplit_option
and pa_HPCdaligner_option. Generally, for large genomes, you should use -s400 (400Mb sequence per block) in
pa_DBsplit_option. This will make a smaller number of jobs but each job will run longer. However, if you have a job
queue system which limits how long a job can run, it might be desirable to have a smaller number for the -s option.

.. _LAsort:

LAsort
++++++

Sort alignment files


FALCON Commands
---------------

.. _fc_run:

fc_run
++++++

This script drives the entire assembly process

.. _fc_consensus:

fc_consensus
++++++++++++

``fc_consensus`` has many options. You can use the parameter :ref:`falcon_sense_option` to control it.
In most cases, the ``--min_cov`` and ``--max_n_read`` are the most important options. ``--min_cov`` controls
when a seed read gets trimmed or broken due to low coverage. ``--max_n_read`` puts a cap on the number of reads
used for error correction. In highly repetitive genome, you will need to make the value for ``--max_n_read``
smaller to make sure the consensus code does not waste time aligning repeats. The longest proper overlaps are used
for correction to reduce the probability of collapsed repeats.

.. _fc_dedup_a_tigs:

fc_dedup_a_tigs
+++++++++++++++
info

.. _fc_graph_to_contig:

fc_graph_to_contig
++++++++++++++++++
info

.. _fc_ovlp_to_graph:

fc_ovlp_to_graph
++++++++++++++++
info

.. _fc_ovlp_filter:

fc_ovlp_filter
++++++++++++++


Troubleshooting FALCON jobs
---------------------------
