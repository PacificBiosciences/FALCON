.. image:: falcon_icon2.png
   :height: 200px
   :width: 200 px
   :alt: Falcon Assembler
   :align: right

.. _tutorial:



Tutorial
========

In this section we will run the FALCON pipeline on an *E. coli* dataset.
We will work through the commands and results and give you ideas of how to assess 
the perfomance of FALCON on your dataset so you can modify parameters and trouble-shoot more 
effectively. This tutorial is a beginners guide to FALCON but assumes bioinformatics fluency.

Input Files
-----------

You will need three types of files to get started, your PacBio data in fasta format (can be one or many files), a 
text file telling FALCON where to find your fasta files, and a :ref:`configuration <Configuration>` file. 
All files except the fasta 
files must be in your job directory.


1. Download *E. coli* dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example *E. coli* dataset can be download from here_. and then unpacked. e.g.:

.. code-block:: bash

	$ wget https://downloads.pacbcloud.com/public/data/git-sym/ecoli.m140913_050931_42139_c100713652400000001823152404301535_s1_p0.subreads.tar.gz
	$ tar -xvzf ecoli.m140913_050931_42139_c100713652400000001823152404301535_s1_p0.subreads.tar.gz 
.. _here: https://downloads.pacbcloud.com/public/data/git-sym/

You should find three fasta files of ~350 Mb each in the newly created subdirectory.


2. Create FOFN
~~~~~~~~~~~~~~

Next, create a "file-of-file-names", ("fofn") with the full path of each fasta file, one per line.

.. code-block:: bash

	/my/path/to/data/ecoli.1.subreads.fasta
	/my/path/to/data/ecoli.2.subreads.fasta
	/my/path/to/data/ecoli.3.subreads.fasta



3. Download configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are running on a cluster with a scheduler use this a starting point: 
:download:`fc_run_ecoli.cfg <cfgs/fc_run_ecoli.cfg>`
If you are running your job locally, try this file: 
:download:`fc_run_ecoli_local.cfg <cfgs/fc_run_ecoli_local.cfg>`

These config files are meant to be starting points only! You will need to make adjustments according
to your particular compute setup.


Note that I have manually specified a seed read length cutoff of 
15kb rather than using an automated cut off (:ref:`length_cutoff <length_cutoff>` = -1, with calculated
length cut off = 22486). The
raw read coverage is very high (>200X); by reducing the seed read length cutoff, we avoid enriching
our seed reads for erroneous chimeric (and very long) reads. Try running the assembly using 
the automated seed read length cut off, you should get a fragmented (28 contigs) and 
incomplete assembly (< 900Mb).

   
Running Falcon
--------------

I send all of my FALCON jobs to the scheduler for ease of tracking job progress. Here is an example
bash script ``run_falcon.sh`` that submits to an SGE_ cluster:

.. code-block:: bash
	
	#!/bin/bash
	#$ -S /bin/bash
	#$ -N myJob
	#$ -cwd
	#$ -q myqueue

	# load dependencies
	module load python/2.7.9 gcc/4.9.2

	# source build
	cd /path/to/my/build/FALCON-integrate
	source env.sh

	# navigate to job directory
	cd /path/to/my/job_dir

	# run it!
	fc_run.py fc_run.cfg


To initiate the FALCON run, I just submit my job to the scheduler with a qsub command:

.. code-block:: bash

	$ qsub run_falcon.sh
	
	
Alternatively, you can add the ``fc_env/bin`` directory to your
``$PATH`` and invoke :doc:`fc_run.py` at the command line with your ``fc_run.cfg`` as the argument.
Note that this shell needs to persist through the entire assembly process so you may need 
to use a window manager like screen_ to maintain your connection.

.. code-block:: bash

    falcon_jobdir$ export PYTHONUSERBASE=/path/to/FALCON-integrate/fc_env
    falcon_jobdir$ export PATH=$PYTHONUSERBASE/bin:$PATH
    falcon_jobdir$ fc_run.py fc_run.cfg


.. _SGE: http://gridscheduler.sourceforge.net/htmlman/manuals.html
.. _screen: https://www.gnu.org/software/screen/manual/screen.html


Assessing Run Progress
----------------------

Refer to the pipeline document for detailed summary of FALCON job directory structure, 
sequence of commands, and output files created.

Counting Completed Jobs
~~~~~~~~~~~~~~~~~~~~~~

The majority of run-time is spent during the daligner phases, performing the alignments and 
then sorting and merging them. To determine how many jobs are performed for each step, refer to ``0-rawreads/run_jobs.sh``.

.. code-block:: bash

    $ grep '^#' run_jobs.sh
    
    	# Daligner jobs (60)
	# Initial sort jobs (400)
	# Check initial .las files jobs (80) (optional but recommended)
	# Remove initial .las files
	# Level 1 merge jobs (20)
	# Check level 2 .las files jobs (20) (optional but recommended)
	# Remove level 1 .las files (optional)

To determine how many jobs have completed, count the sentinel files that indicate a job is complete.
For example:

.. code-block:: bash

	$ find 0-rawreads/ -name "job*done" | wc -l
	60
	
	$ find 0-rawreads/ -name "m_*done" | wc -l
	20


Assessing Run Performance
-------------------------

Raw and Pread Coverage and Quality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The *E. coli* subreads
are a total of 1.01 Gb of data in 105,451 reads. (:download:`countFasta.pl <countFasta.pl>` 
is a useful script by Joseph Fass and Brad Sickler at UC Davis for calculating total sequence
length and other assembly metrics).

You can confirm that your dazzler database was correctly constructed using a utility from the dazzler_ suite:

.. _dazzler: https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide/

.. code-block:: bash 

	$ DBstats raw_reads.db > raw_reads.stats
	$ head raw_reads.stats
	
	Statistics for all reads of length 500 bases or more
	
		90,747 reads        out of         105,451  ( 86.1%)
		964,281,429 base pairs   out of   1,013,118,375  ( 95.2%)
	
		10,626 average read length
		6,805 standard deviation
	
	Base composition: 0.248(A) 0.242(C) 0.263(G) 0.246(T)
	
	Distribution of Read Lengths (Bin size = 1,000)
	
		Bin:      Count  % Reads  % Bases     Average
		45,000:     1      0.0      0.0       45611


You can see that we discarded 13.9% of the raw bases and 4.8% of the reads by employing a 
raw read length cut off of 500bp in the :ref:`DBsplit <dazzDBsplit>` options. This file can
also be used to plot a :ref:`histogram <RHists>` of raw read lengths.

The genome of this *E. coli* strain is 4.65 Mb long for a raw read coverage of ~207 fold.
Confirm this with the preassembly report:

.. code-block:: bash 

	$ cat 0-rawreads/report/pre_assembly_stats.json
	
	"genome_length": 4652500,
	"length_cutoff": 15000,
	"preassembled_bases": 350302363, 	
	"preassembled_coverage": 75.293,	
	"preassembled_mean": 10730.33,		
	"preassembled_n50": 16120,			
	"preassembled_p95": 22741,
	"preassembled_reads": 32646,
	"preassembled_seed_fragmentation": 1.451,	# number split preads / seed reads 
	"preassembled_seed_truncation": 4453.782,	# ave bp lost per pread due to low cov
	"preassembled_yield": 0.758,			# total pread bp / seed read bp
	"raw_bases": 964281429,
	"raw_coverage": 207.261,
	"raw_mean": 10626.042,
	"raw_n50": 14591,
	"raw_p95": 23194,
	"raw_reads": 90747,
	"seed_bases": 461851093,	
	"seed_coverage": 99.269,			# raw base coverage depth on seed reads
	"seed_mean": 20029.103,
	"seed_n50": 19855,
	"seed_p95": 28307,
	"seed_reads": 23059

A note on these statistics: in the process of created preads, seeds reads with insufficient
raw read coverage (usually due to base errors) will be split or truncated. The preassembled seed
fragmentation, truncation, and yield stats summarize the quality of pread assembly. 
A good preassembled yield should be greater than 50%. Note that if an automated seed read length
is used for this data (22486), preassembled seed read truncation is ~6kb, indicating that many of the longest
raw reads are not supported by the rest of the data.

You can similarly summarize the contents of the dazzler database for preads using DBstats 
and plotting in :ref:`R <RHists>`.

Contig Stats
~~~~~~~~~~~~

When your run is complete, you can summarize your assembly stats using the countFasta.pl script:

 .. code-block:: bash
	
	$ countFasta.pl p_ctg.fa > p_ctg.stats
	$ countFasta.pl a_ctg.fa > a_ctg.stats
	$ tail p_ctg.stats
	
	Total length of sequence:	4635395 bp
	Total number of sequences:	1
	N25 stats:			25% of total sequence length is contained in the 1 sequences >= 4635395 bp
	N50 stats:			50% of total sequence length is contained in the 1 sequences >= 4635395 bp
	N75 stats:			75% of total sequence length is contained in the 1 sequences >= 4635395 bp
	Total GC count:			2352187 bp
	GC %:				50.74 %


Assembly Graph and Pread Overlaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assembly contiguity can be enhanced by adjusting a few parameters in the last stage of the 
assembly process. You can try a grid of :ref:`pread length cut offs <length_cutoff_pr>` for 
the filtering of the final
overlaps in the assembly graph. In a general sense, longer pread length cut offs will increase the
 contiguity (contig N50) in your assembly, but may result in shorter over all assembly length. 
To try different length cut off, rename your 2-asm-falcon dir,
modify the config file, rename the log and mypwatcher directory, and restart falcon:

.. code-block:: bash
	
	$ mv 2-asm-falcon 2-asm-falcon_12kb
	$ mv mypwatcher/ mypwatcher0/
	$ mv all.log all0.log
	$ qsub run_falcon.sh


The other parameter to adjust is the number of overlaps in the assembly graph. First, look
at a histogram of the number of overlaps on the 5' and 3' end of each read. Run the falcon utility:

.. code-block:: bash
	
	$ cd 2-asm-falcon
	$ fc_ovlp_stats.py --fofn 1-preads_ovl/merge-gather/las.fofn > ovlp.stats
	
Then plot :ref:`histograms <OvlpHists>` of the number of 5' and 3' overlaps between preads in R.
This can inform your parameters for :ref:`sge_option_fc <sge_option_fc>` where ``min_cov`` and ``max_cov``
should flank the bulk of the distribution. For repetative genomes, a second mode in the :ref:`distribution <RepeatOvlp>`
may appear, containing preads ending or begining in repetative material. It is best to choose a ``max_cov``
to the left of the repeat mode that removes these repetative overlaps.



Troubleshooting Run
-------------------

If you find your run has died here are some suggestions of how to restart,
in order of increasing difficulty:

Simple Restart
~~~~~~~~~~~~~~

Simply rename your log file and ``mypwatcher`` directory and restart the pipeline. Renaming these
files preserves them for you reference, and by removing the original mypwatcher directory
the pipeline, when restarted, will scan your job directory for completed jobs and pick up where it left off:

.. code-block:: bash

	$ mv mypwatcher/ mypwatcher0/
	$ mv all.log all0.log
	$ qsub run_falcon.sh


Directory Cleanup and Restart
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, determine which job caused the run to fail. For example:

.. code-block:: bash
	
	$ grep 'ERROR' all.log
	
	2016-11-21 03:21:39,482 - pypeflow.simple_pwatcher_bridge - ERROR - Task Node(0-rawreads/m_00210) failed with exit-code=99
	2016-11-21 03:21:39,482 - pypeflow.simple_pwatcher_bridge - ERROR - Failed to clean-up FakeThread: jobid=Pcfbdb8b3c50d5e status='EXIT '

Delete all directories that failed, then rename the log file and ``mypwatcher`` as above:

.. code-block:: bash

	$ rm -rf 0-rawreads/m_00210
	$ mv mypwatcher/ mypwatcher0/
	$ mv all.log all0.log
	$ qsub run_falcon.sh

You can find out more details about the failed jobs in ``mypwatcher/`` to diagnose the problem.

.. code-block:: bash

	$ less mypwatcher/jobs/Pcfbdb8b3c50d5e/stderr
	$ less mypwatcher/jobs/Pcfbdb8b3c50d5e/stdout


Manual Running of Failed Jobs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your job still fails, try manually running the problematic jobs. Search in the job
directory for the shell script containing the individual tasks and try manually running 
the shell script or individual tasks:

.. code-block:: bash

	$ ls job_0000
	
	job_0000_done  L1.19.5.las  L1.19.7.las  L1.5.19.las  L1.7.19.las   raw_reads.db  run.sh       task.json
	L1.19.4.las    L1.19.6.las  L1.4.19.las  L1.6.19.las  pwatcher.dir  rj_0000.sh    run.sh.done  task.sh

	$ head job_0000/rj_0000.sh -n 12

	#!/bin/bash
	set -vex

	db_dir=/lustre/hpcprod/skingan/FALCON_tutorial/ecoli/0-rawreads
	ln -sf ${db_dir}/.raw_reads.bps .
	ln -sf ${db_dir}/.raw_reads.idx .
	ln -sf ${db_dir}/raw_reads.db .
	ln -sf ${db_dir}/.raw_reads.dust.anno .
	ln -sf ${db_dir}/.raw_reads.dust.data .
	daligner -v -t16 -H22486 -e0.7 -s1000 raw_reads.19 raw_reads.4 raw_reads.5 raw_reads.6 raw_reads.7
	LAcheck -v raw_reads *.las
	LAsort -v raw_reads.4.raw_reads.19.C0 raw_reads.4.raw_reads.19.N0 raw_reads.4.raw_reads.19.C1 raw_reads.4.raw_reads.19.N1 raw_reads.4.raw_reads.19.C2 raw_reads.4.raw_reads.19.N2 raw_reads.4.raw_reads.19.C3 raw_reads.4.raw_reads.19.N3 && LAmerge -v L1.4.19 raw_reads.4.raw_reads.19.C0.S raw_reads.4.raw_reads.19.N0.S raw_reads.4.raw_reads.19.C1.S raw_reads.4.raw_reads.19.N1.S raw_reads.4.raw_reads.19.C2.S raw_reads.4.raw_reads.19.N2.S raw_reads.4.raw_reads.19.C3.S raw_reads.4.raw_reads.19.N3.S

Once these jobs have run to completion, you can try restarting the pipeline.
