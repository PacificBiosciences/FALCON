.. _pipeline:

.. image:: falcon_icon2.png
   :height: 200px
   :width: 200 px
   :alt: Falcon Assembler
   :align: right


########
Pipeline
########

FALCON
======

A FALCON job can be broken down into 3 steps:

1. :ref:`Overlap detection and error correction of rawreads <step_one>`
2. :ref:`Overlap detection between corrected reads <step_two>`
3. :ref:`String Graph assembly of corrected reads <step_three>`

Each step is performed in it's own subdirectory within the FALCON job

.. code-block:: bash

    falcon_job/
        ├── 0-rawreads     # Raw read error correction directory
        ├── 1-preads_ovl   # Corrected read overlap detection
        ├── 2-asm-falcon   # String Graph Assembly
        ├── scripts
        └── sge_log        # Job scheduler logs

The assembly process is driven by the script :ref:`fc_run.py` which should be sent to the scheduler or run on a head node
as it needs to persist throughout the entire assembly process.
It takes as input a single :ref:`config <Configuration>` file typically named :ref:`fc_run.cfg`, which references a list
of fasta input files. The config file can be configured to
run locally, or submit to a job scheduler. However, if your dataset is anything larger than a bacterial sized
genome and unless you've tuned your system specifically for the organism you're trying to assemble, then most likely you
should be running on a cluster in order to more effectively leverage your computational resources.

The configuration file also allows you to control other aspects of your job such as how your compute resources are
distributed as well as set many parameters to help you reach an "optimized" assembly according to the nature of
your input data. Unfortunately at this point there is no "magic" way to auto-tune the parameters so you should
probably spend some time in the :ref:`configuration` section to understand what options are available to you.
Some example configuration files can be found :ref:`here <configuration>`


.. _step_one:

Step 1: Overlap detection and error correction of raw reads
-----------------------------------------------------------

The first step of the pipeline is to identify all overlaps in the raw reads. Currently this is performed with
a modified version of Gene Myers' DALIGNER_.

In order to identify overlaps, your :term:`raw reads` must first be converted from fasta format into a dazzler
database. This is a very I/O intensive process and will be run from the node where ``fc_run.py`` was executed. If this
is an issue, you should submit the command with a wrapper script to your grid directly.

Once the database has been created and partitioned according to the parameters set in your
:ref:`fc_run.cfg <configuration>`, an all vs all comparison of the reads must be performed. Accordingly, due to the
all vs all nature of the search this is the most time consuming step in the assembly process. To walk through the
actual steps of this part of the pipeline you should take a look at `0-rawreads/prepare_rdb.sub.sh`.
Essentially it consists of running:

1. :ref:`fasta2DB <fasta2DB>` to format the database
2. :ref:`DBsplit <DBsplit>` to partition the database
3. :ref:`HPC.daligner <HPC.daligner>` to generate the :ref:`daligner` commands necessary for all-vs-all comparison

After overlaps have been detected, you will be left with many ``job_*`` directories full of alignment files ``*.las`` 
containing the information about the overlaps. After merging the alignment files (see ``m_*`` directories), the 
next step is to error correct the reads leveraging the overlap information. In the ``0-rawreads/preads`` directory you 
will find a series of scripts for
performing the error correction. The process basically consists of using ``LA4Falcon`` with a length cutoff and piping the
output to :ref:`fc_consensus.py <fc_consensus>` to generate a fasta file with corrected reads.


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

The following parameters affect this step directly:

* :ref:`sge_option_da <sge_option_da>`
* :ref:`sge_option_la <sge_option_la>`
* :ref:`pa_concurrent_jobs <pa_concurrent_jobs>`
* :ref:`cns_concurrent_jobs <cns_concurrent_jobs>`
* :ref:`pa_DBsplit_option <pa_DBsplit_option>`
* :ref:`falcon_sense_option <falcon_sense_option>`

.. _DALIGNER: http://dazzlerblog.wordpress.com
.. _Dazzler: https://dazzlerblog.wordpress.com/2014/06/01/the-dazzler-db/


.. _step_two:

Step 2: Overlap detection of corrected reads
--------------------------------------------

The only conceptual difference between the first and second overlap steps is that consensus calling is
not performed in the second step. After :term:`pread` overlap detection, it's simply a
matter of extracting the information from the corrected reads database with ``DB2Falcon -U preads``.

Depending on how well the error-correction step proceeded as well as the how much
initial coverage was fed into the pipeline (e.g. :ref:`length_cutoff <length_cutoff>`), the input data for this 
step should be significantly reduced and thus, the second overlap detection step 
will proceed significantly faster.

The commands in this step of the pipeline are very similar to before albeit with different parameter settings to account
for the reduced error-rate of the :term:`preads <pread>`. See the driver script ``prepare_pdb.sub.sh`` for 
details on actual parameter settings used.

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

The following parameters affect this step directly:

* :ref:`sge_option_pda <sge_option_pda>`
* :ref:`sge_option_pla <sge_option_pla>`
* :ref:`ovlp_concurrent_jobs <ovlp_concurrent_jobs>`
* :ref:`ovlp_DBsplit_option <ovlp_DBsplit_option>`
* :ref:`ovlp_HPCdaligner_option <ovlp_HPCdaligner_option>`


.. _step_three:

Step 3: String Graph assembly
-----------------------------

The final step of the FALCON Assembly pipeline is generation of the final :term:`String Graph` assembly and 
output of contig sequences in
fasta format. Four commands are run in the final phase of FALCON:

1. :ref:`fc_ovlp_filter <fc_ovlp_filter>` Filters overlaps based on the criteria provided in :ref:`fc_run.cfg`
2. :ref:`fc_ovlp_to_graph <fc_ovlp_to_graph>` constructs an overlap graph of reads larger than the length cutoff
3. :ref:`fc_graph_to_contig <fc_graph_to_contig>` generates fasta files for contigs from the overlap graph.
4. :ref:`fc_dedup_a_tigs <fc_dedup_a_tigs>` removes duplicate associated contigs

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

The following parameters affect this step directly:

* :ref:`sge_option_fc <sge_option_fc>`
* :ref:`overlap_filtering_setting <overlap_filtering_setting>`
* :ref:`length_cutoff_pr <length_cutoff_pr>`


FALCON_unzip
============

FALCON_unzip operates from a completed FALCON job directory. There are x steps to the FALCON_unzip pipeline

1. Read tracking