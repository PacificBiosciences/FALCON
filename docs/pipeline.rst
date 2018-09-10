.. _pipeline:

.. caution:: These documents refer to an obsolete way of installing and running FALCON. They will remain up for historical context and for individuals still using the older version of FALCON/FALCON_unzip.

.. attention:: The current PacBio Assembly suite documentation which includes new bioconda instructions for installing FALCON, FALCON_unzip and their associated dependencies can be found here `pb_assembly <http://github.com/PacificBiosciences/pb-assembly>`_

.. image:: media/falcon_icon2.png
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
        ├── 0-rawreads/     # Raw read error correction directory
        ├── 1-preads_ovl/   # Corrected read overlap detection
        ├── 2-asm-falcon/   # String Graph Assembly
        ├── mypwatcher/     # Job scheduler logs
        ├── scripts/
        └── sge_log/        # deprecated

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
actual steps of this part of the pipeline you should take a look at ``0-rawreads/prepare_rdb.sub.sh``.
Essentially it consists of running:

1. :ref:`fasta2DB <dazzfasta2DB>` to format the database
2. :ref:`DBsplit <dazzDBsplit>` to partition the database
3. :ref:`HPC.daligner <dazzHPC.daligner>` to generate the :ref:`daligner` commands necessary for all-vs-all comparison

After overlaps have been detected, you will be left with many ``job_*`` directories full of alignment files ``*.las`` 
containing the information about the overlaps. After merging the alignment files (see ``m_*`` directories), the 
next step is to error correct the reads leveraging the overlap information. In the ``0-rawreads/preads`` directory you 
will find a series of scripts for
performing the error correction. The process basically consists of using ``LA4Falcon`` with a length cutoff and piping the
output to :ref:`fc_consensus.py <fc_consensus>` to generate a fasta file with corrected reads.


.. code-block:: bash

    0-rawreads/
        ├── job_*                     # dirs for all of the daligner jobs
        ├── m_*/                      # dirs for all of the LA4Merge jobs
        ├── preads/                   # sub-dir for preads generation
        ├── report/		              # pre-assembly stats
        ├── cns-scatter/	          # dir of scripts for falcon-consensus jobs
        ├── daligner-scatter/	      # dir of scripts for daligner jobs
        ├── merge-scatter/	          # dir of scripts for LAMerge jobs
        ├── merge-gather/	          # dir of scripts for gathering LAMerge inputs
        ├── raw-gather/	      	      # dir of scripts for gathering daligner jobs for merging
        ├── input.fofn               # list if your input *.fasta files
        ├── length_cutoff             # text file with length cutoff for seed reads
        ├── pwatcher.dir	          # dir of individual pipeline jobs stderr and stdout
        ├── prepare_rdb.sh            # env wrapper script
        ├── raw_reads.db              # dazzler DB file
        ├── raw-fofn-abs	          # dir of scripts for gathering raw reads inputs
        ├── rdb_build_done            # database construction sentinel file
        ├── run_jobs.sh              # listing of all overlap step commands
        ├── run.sh		              # masker job script
        ├── run.sh.done		          # sentinel file for all jobs
        ├── task.json		         # json file specifying inputs, outputs, and params
        └── task.sh		             # script to run json file



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
        ├── job_*/ 		    # directories for daligner jobs
        ├── m_*/                    # directories for LA4Merge jobs
        ├── db2falcon/	      	    # dir of scripts for formatting preads for falcon
        ├── gathered-las/	    # dir of scripts for gathering daligner jobs
        ├── merge-gather/	    # dir of scripts for gathering LAMerge inputs
        ├── merge-scatter/	    # dir of scripts for LAMerge jobs
        ├── daligner-scatter/	    # dir of scripts for daligner jobs
        ├── pdb_build_done          # sentinel file for pread DB building
        ├── preads.db               # preads dazzler DB
        ├── prepare_pdb.sh          # env wrapper script
        ├── pwatcher.dir	    # dir of individual pipeline jobs stderr and stdout
        ├── run_jobs.sh             # listing of all pread overlap job commands
        ├── run.sh		    # masker job script
        ├── run.sh.done		    # sentinel file for all jobs
        ├── task.json		    # json file specifying inputs, outputs, and params
        └── task.sh		    # script to run json file

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
output of contig sequences in fasta format. Four commands are run in the final phase of FALCON:

1. :ref:`fc_ovlp_filter <fc_ovlp_filter.py>` - Filters overlaps based on the criteria provided in :ref:`fc_run.cfg`
2. :ref:`fc_ovlp_to_graph <fc_ovlp_to_graph.py>` - Constructs an overlap graph of reads larger than the length cutoff
3. :ref:`fc_graph_to_contig <fc_graph_to_contig.py>` - Generates fasta files for contigs from the overlap graph.
4. :ref:`fc_dedup_a_tigs <fc_dedup_a_tigs.py>` - Removes duplicate associated contigs

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
        ├── run_falcon_asm.sh            # env wrapper script
        ├── task.json		              # json file specifying inputs, outputs, and params
        ├── task.sh		                  # script to run json file
        ├── run.sh.done		              # sentinel file for all jobs
        └── run.sh                       # Assembly driver script


The following parameters affect this step directly:

* :ref:`sge_option_fc <sge_option_fc>`
* :ref:`overlap_filtering_setting <overlap_filtering_setting>`
* :ref:`length_cutoff_pr <length_cutoff_pr>`


FALCON_unzip
============

`FALCON_unzip`_ operates from a completed FALCON job directory. After tracking the raw reads to contig,
A FALCON_unzip job can be broken down into 3 steps

1. :ref:`Identify SNPs and assign phases <unzip_step_one>`
2. :ref:`Annotate Assembly graph with Phases <unzip_step_two>`
3. :ref:`Graph building <unzip_step_three>`

.. code-block:: bash

        3-unzip/
        ├── 0-phasing/                  # Contig phasing jobs
        ├── 1-hasm/                     # Contig Graph assembly information
        ├── read_maps/                  # rawread_to_contigs; read_to_contig_map
        ├── reads/                      # raw read fastas for each contig
        ├── all_p_ctg.fa                # partially phased primary contigs
        ├── all_h_ctg.fa                # phased haplotigs
        ├── all_p_ctg_edges             # primary contig edge list
        ├── all_h_ctg_edges             # haplotig edge list
        ├── all_h_ctg_ids               # haplotig id index
        └── all_phased_reads            # table of all phased raw reads


.. _FALCON_unzip:: https://github.com/PacificBiosciences/FALCON_unzip

.. _unzip_step_one:

Step 1: Identify SNPs and assign phases
---------------------------------------

Inside of ``0-phasing/`` you vill find a number of directories for each contig. Each contains the scripts
to map the raw reads to the contigs and subsequently identify SNPs. The generated SNP tables can
subsequently be used to assign phases to reads.


.. _unzip_step_two:

Step 2: Graph annotation and haplotig
-------------------------------------

Inside of ``1-hasm/`` you can find the driver script ``hasm.sh`` which contains the commands necessary to
filter overlaps and traverse the assembly graph paths and subsequently output phased contig sequence.
Assembly Graphs for each contig as well as fasta files for the partially phased primary contigs and fully phased
haplotigs can be found in each ``1-hasm/XXXXXXF`` directory.


.. _unzip_step_three:

Step 3: Call Consensus (Optional)
---------------------------------

Finally, the ``FALCON_unzip`` pipeline can optionally be used to run quiver and call high quality consensus. This step
takes as input the primary contig and haplotig sequences output in the previous step. For convenience, these files
have all been concatenated together into ``3-unzip/all_p_ctg.fa`` and ``3-unzip/all_h_ctg.fa`` respectively.
The final consensus output can be found in ``falcon_jobdir/4-quiver/cns_output/*.fast[a|q]``.
In order to run the consensus step as part of the FALCON_unzip pipeline, You need to provide the :ref:`input_bam_fofn`
:ref:`fc_unzip.cfg` option in order for this to work.

