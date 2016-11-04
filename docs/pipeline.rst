.. _pipeline:

FALCON Pipeline
===============


A FALCON job can be broken down into 3 steps:

1. Overlap detection and :term:`error correction` of rawreads
2. Overlap detection between corrected reads (:term:`pread`)
3. :term:`String Graph` assembly of corrected reads

Each step is performed in it's own subdirectory within the FALCON job

.. code-block:: bash

    falcon_job/
        ├── 0-rawreads
        ├── 1-preads_ovl
        ├── 2-asm-falcon
        ├── scripts
        └── sge_log

The assembly process is driven by the script ``fc_run.py`` which should be run an a head node as it needs to persist
throughout the entire assembly process.
It takes as input a single config file typically named `fc_run.cfg` and can be configured to run locally, or submit
to a job scheduler



Overlap detection and error correction of raw reads (0-rawreads)
----------------------------------------------------------------

The first step of the pipeline is to identify all overlaps in the raw reads. Currently this is all performed with
a modified version of Gene Myers' Daligner_.

In order to identify overlaps, your :term:`raw reads` must first be converted from fasta format into a Dazzler_ database.
Once the database has been created and partitioned according to the parameters set in your ``fc_run.cfg``,
an all vs all comparison of the reads must be performed. Accordingly, due to the all vs all nature of the search
this is the most time consuming step in the assembly process. To walk through the actual steps of this part of the
pipeline you should take a look at `0-rawreads/prepare_rdb.sub.sh`. Essentially it consists of running:

1. ``fasta2DB`` to format the database
2. ``DBSplit`` to partition the database
3. ``HPC.daligner`` to generate the ``daligner`` commands necessary for all-vs-all comparison

After overlaps have been detected, you will be left with many directories full of alignment files ``*.las`` containing
the information about the overlaps. After merging the alignment files, the next step is to error correct the reads
leveraging the overlap information. In the ``0-rawreads/preads`` directory you will find a series of scripts for
performing the error correction. The process basically consists of using ``LA4Falcon`` with a length cutoff and piping the
output to ``fc_consensus`` to generate a fasta file with corrected reads. ``fc_consensus.py`` has many options.
You can use the parameter :ref:`falcon_sense_option` to control it. In most cases, the ``--min_cov`` and ``--max_n_read``
are the most important options. ``--min_cov`` controls when a seed read gets trimmed or broken due to low coverage.
``--max_n_read`` puts a cap on the number of reads used for error correction. In highly repetitive genome, you will
need to make the value for ``--max_n_read`` smaller to make sure the consensus code does not waste time aligning
repeats. The longest proper overlaps are used for correction to reduce the probability of collapsed repeats.

One can use :ref:`cns_concurrent_jobs` to control the maximum number of concurrent jobs submitted to the job management
system. The ``out.XXXXX.fasta`` files produced are used as input for the next step in the pipeline.

.. code-block:: bash

    0-rawreads/
        ├── pre_assembly_stats.json  #Pre-assembly stats
        ├── cns_done                  #consensus sentinel file
        ├── preads/                   #sub-dir for error correction
        ├── m_*/                      #All of the merge jobs
        ├── da_done                   # daligner sentinel file
        ├── job_0*                   #All of the daligner jobs
        ├── length_cutoff            #text file with just the length cutoff
        ├── raw_reads.db             #Dazzler DB file
        ├── rdb_build_done           #database construction sentinel file
        ├── run_jobs.sh              #Listing of all overlap step commands
        ├── input.fofn               #List if your input *.fasta files
        ├── prepare_rdb.sh           #env wrapper script
        └── prepare_rdb.sub.sh       #Driver script for this step in the pipeline



.. _Daligner: http://dazzlerblog.wordpress.com
.. _Dazzler: https://dazzlerblog.wordpress.com/2014/06/01/the-dazzler-db/


Overlap detection of corrected reads (1-preads_ovl)
---------------------------------------------------

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
        ├── db2falcon_done
        ├── preads4falcon.fasta
        ├── run_db2falcon.sh
        ├── run_db2falcon.sub.sh
        ├── p_merge_done
        ├── m_*/
        ├── da_done
        ├── job_*/
        ├── pdb_build_done
        ├── preads.db¡
        ├── run_jobs.sh
        ├── prepare_pdb.sh
        ├── prepare_pdb.sub.sh
        └── input_preads.fofn


Graph assembly (2-asm_falcon)
-----------------------------

The final step of the FALCON Assembly pipeline is generation of the final graph assembly and output in fasta format.
There are 4 commands being run in the final phase of the FALCON assembly pipeline:

1. ``fc_ovlp_filter`` is responsible for filtering your overlaps based on the criteria you provided in fc_run.cfg
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
        ├── ctg_paths                    # Read paths for each contig
        ├── fc_ovlp_to_graph.log        # logfile for process of converting overlaps to assembly graph
        ├── utg_data                     #
        ├── sg_edges_list               # list of all edges
        ├── chimers_nodes               #
        ├── preads.ovl                  # List of all overlaps between preads
        ├── las.fofn                    # List of *.las files for input
        ├── run_falcon_asm.sh           # env wrapper script
        └── run_falcon_asm.sub.sh       # Assembly driver script



Troubleshooting FALCON jobs
---------------------------
