.. image:: falcon_icon2.png
   :height: 200px
   :width: 200 px
   :alt: Falcon Assembler
   :align: right


.. _parameters:

##########
Parameters
##########


.. _configuration:

Configuration
=============

Here are some example ``fc_run.cfg`` and ``fc_unzip.cfg`` files. We make no guarantee that they will work with your
dataset and cluster configuration. We merely provide them as starting points that have proven themselves on internal
datasets. A lot of your success will depend purely on the quality of the input data prior to even engaging the FALCON
pipeline. Also, these particular configs were designed to work in our SGE compute cluster, so some tuning will likely
be necessary on your part. You should consult with your HPC administrator to assist in tuning to your cluster.

FALCON Parameter sets
---------------------

:download:`fc_run_fungal.cfg <cfgs/fc_run_fungal.cfg>` - Has worked well on a 40Mb fungal genome

:download:`fc_run_human.cfg <cfgs/fc_run_human.cfg>` - Has worked well on at least one human dataset

:download:`fc_run_bird.cfg <cfgs/fc_run_bird.cfg>` - Has worked well on at least one avian dataset

:download:`fc_run_yeast.cfg <cfgs/fc_run_yeast.cfg>` - Has worked well on at least one yeast dataset

:download:`fc_run_dipteran.cfg <cfgs/fc_run_dipertan.cfg>` - Has worked well on at least one dipteran (insect) dataset

:download:`fc_run_mammal.cfg <cfgs/fc_run_mammal.cfg>` - Has worked well on at least one mammalian dataset

:download:`fc_run_plant.cfg <cfgs/fc_run_plant.cfg>` - Has worked well on at least one plant (Ranunculales) dataset

:download:`fc_run_arabidopsis.cfg <cfgs/fc_run_arabidopsis.cfg>` - Configuration for arabidopsis assembly in Chin et al. 
2016 


FALCON_unzip Parameter sets
---------------------------

:download:`fc_unzip.cfg <cfgs/fc_unzip.cfg>` - General all purpose unzip config


Available Parameters
====================

.. _fc_run.cfg:

fc_run.cfg
----------

.. _input_fofn:

input_fofn <str>
   filename for the file-of-filenames (fofn)
   Each line is fasta filename.
   Any relative paths are relative to the location of the input_fofn.

.. _input_type:

input_type <str>
   "raw" or "preads"


.. _genome_size:

genome_size <int>
   estimated number of base-pairs in haplotype

.. _seed_coverage:

seed-coverage <int>
   requested coverage for auto-calculated cutoff

.. _length_cutoff:

length_cutoff <int>
   minimum length of seed-reads used for pre-assembly stage
   If '-1', then auto-calculate the cutoff based on genome_size and seed_coverage.

.. _length_cutoff_pr:

length_cutoff_pr <int>
   minimum length of seed-reads used after pre-assembly, for the "overlap" stage


.. _target:

target <str>
   "assembly" or "preads"
   If "preads", then pre-assembly stage is skipped and input is assumed to be preads.


.. _default_concurrent_jobs:

default_concurrent_jobs <int>
   maximum concurrency
   This applies even to "local" (non-distributed) jobs.

.. _pa_concurrent_jobs:

pa_concurrent_jobs <str>
   Concurrency settings for pre-assembly

.. _cns_concurrent_jobs:

cns_concurrent_jobs <str>
   Concurrency settings for consensus calling

   One can use cns_concurrent_jobs to control the maximum number of concurrent consensus jobs submitted to the
   job management system. The ``out.XXXXX.fasta`` files produced are used as input for the next step in the pipeline.


.. _ovlp_concurrent_jobs:

ovlp_concurrent_jobs <str>
   Concurrency settings for Overlap detection

.. _job_type:

job_type <str>
   grid submission system, or "local"
   Supported types include: "sge", "lsf", "pbs", "torque", "slurm", "local"
   case-insensitive

.. _job_queue:

job_queue <str>
   grid job-queue name
   Can be overridden with section-specific sge_option_*

.. _sge_option_da:

sge_option_da <str>
   Grid concurrency settings for initial daligner steps ``0-rawreads/``

.. _sge_option_la:

sge_option_la <str>
   Grid concurrency settings for initial las-merging ``0-rawreads/``

.. _sge_option_cns:

sge_option_cns <str>
   Grid concurrency settings for error correction consensus calling

.. _sge_option_pda:

sge_option_pda <str>
   Grid concurrency settings for daligner on preads ``1-preads_ovl/``

.. _sge_option_pla:

sge_option_pla <str>
   Grid concurrency settings for las-merging on preads in ``1-preads_ovl/``

.. _sge_option_fc:

sge_option_fc <str>
   Grid concurrency settings for stage 2 in ``2-asm-falcon/``

.. _pa_DBdust_option:

pa_DBdust_option <str>
   Passed to ``DBdust``. Used only if ``dust = true``.

.. _pa_DBsplit_option:

pa_DBsplit_option <str>
   Passed to ``DBsplit`` during pre-assembly stage.


.. _pa_HPCdaligner_option:

pa_HPCdaligner_option <str>
   Passed to ``HPC.daligner`` during pre-assembly stage.
   We will add ``-H`` based on``length_cutoff``.

   The ``-dal`` option also controls the number of jobs being spawned. The number
   for the ``-dal`` option determines how many blocks are compared to each in single jobs. Having a larger number
   will spawn a fewer number of larger jobs, while the opposite will give you a larger number of small jobs. This
   will depend on your on your compute resources available.

   In this workflow, the trace point generated by ``daligner`` is not used. ( Well, to be efficient, one should use the trace
   points but one have to know how to pull them out correctly first. ) The ``-s1000`` argument makes the trace points
   sparse to save some disk space (not much though). We can also ignore all reads below a certain
   threshold by specifying a length cutoff with ``-l1000``.

   The biggest difference between this parameter and the ``ovlp_HPCdaligner_option`` parameter is that the latter needs
   to have a relaxed error rate switch ``-e`` as the alignment is being performed on uncorrected reads.

.. _pa_dazcon_option:

pa_dazcon_option <str>
   Passed to ``dazcon``. Used only if ``dazcon = true``.

.. _falcon_sense_option:

falcon_sense_option <str>
   Passed to ``fc_consensus``.
   Ignored if ``dazcon = true``.

.. _falcon_sense_skip_contained:

falcon_sense_skip_contained <str>
   Causes ``-s`` to be passed to ``LA4Falcon``. Rarely needed.

.. _ovlp_DBsplit_option:

ovlp_DBsplit_option <str>
   Passed to ``DBsplit`` during overlap stage.

.. _ovlp_HPCdaligner_option:

ovlp_HPCdaligner_option <str>
   Passed to ``HPC.daligner`` during overlap stage.

.. _overlap_filtering_setting:

overlap_filtering_setting <str>
   Passed to ``fc_ovlp_filter`` during assembly stage.

.. _fc_ovlp_to_graph_option:

fc_ovlp_to_graph_option <str>
   Passed to ``fc_ovlp_to_graph``.

.. _skip_checks:

skip_check <bool>
   If "true", then skip ``LAcheck`` during ``LAmerge``/``LAsort``.
   (Actually, ``LAcheck`` is run, but failures are ignored.)
   When ``daligner`` bugs are finally fixed, this will be unnecessary.


.. _dust:

dust <bool>
   If true, then run ``DBdust`` before pre-assembly.

.. _dazcon:

dazcon <bool>
   If true, then use ``dazcon`` (from pbdagcon repo).


.. _stop_all_jobs_on_failure:

stop_all_jobs_on_failure <bool>
   DEPRECATED
   This was used for the old pypeFLOW refresh-loop, used by ``run0.py``.
   (This is *not* the option to let jobs currently in SGE (etc) to keep running, which is still TODO.)

.. _use_tmpdir:

use_tmpdir <bool>
   (boolean string) whether to run each job in ``TMPDIR`` and copy results back to nfs
   If "true", use ``TMPDIR``. (Actually, ``tempfile.tmpdir``. See standard Python docs: https://docs.python.org/2/library/tempfile.html )
   If the value looks like a path, then it is used instead of ``TMPDIR``.


.. _fc_unzip.cfg:

fc_unzip.cfg
------------

job_type <str>
   same as above.
   grid submission system, or "local"
   Supported types include: "sge", "lsf", "pbs", "torque", "slurm", "local"
   case-insensitive

input_fofn <str>
   This will be the same input file you used in your :ref:`fc_run.cfg`

.. _input_bam_fofn:

input_bam_fofn <str>
   List of movie bam files. Only necessary if performing consensus calling step at the end.

smrt_bin <str>
   path to ``bin`` directory containing samtools, blasr, and various GenomicConsensus utilities

jobqueue <str>
   Queue to submit SGE jobs to.

sge_phasing <str>
   Phasing grid settings. Example: ``-pe smp 12 -q %(jobqueue)s``

sge_quiver <str>
   Consensus calling grid settings. Example ``-pe smp 24 -q %(jobqueue)s``

sge_track_reads <str>
   Read tracking grid settings. Example ``-pe smp 12 -q %(jobqueue)s``

sge_blasr_aln <str>
   ``blasr`` alignment grid settings. Example ``-pe smp 24 -q %(jobqueue)s``

sge_hasm <str>
   Final haplotyped assemble grid settings Example ``-pe smp 48 -q %(jobqueue)s``

unzip_concurrent_jobs <int>
   Number of concurrent unzip jobs to run at a time

quiver_concurrent_jobs <int>
   Number of concurrent consensus calling jobs to run