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
the perfomance of FALCON on your dataset so you can modify parameters and trouble-shoot more effectively.

Input Files
-----------

You will need three types of files to get started, your PacBio data in fasta format (can be one or many files), a 
text file telling FALCON where to find your fasta files, and a configuration file. All files except the fasta 
files must be placed in your job directory. More about the structure of FALCON job directories can be found in the 
:ref:`pipeline` section.


1. Download *E. coli* dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The *E. coli* example dataset is provided with the FALCON-integrate distribution. To download it, first
follow the :ref:`quick_start` to install FALCON and set up your environment. If you have already installed 
FALCON, you just need to load any dependencies, then run ``source env.sh`` command from inside the 
FALCON-integrate directory. Next, navigate to the FALCON-integrate/FALCON-examples directory, and download the dataset, 
e.g.:

.. code-block:: bash

	module load python/2.7.9 gcc/4.9.2 git 
	cd FALCON-integrate
	source env.sh
	cd FALCON-example
	make setup-ecoli
	
You should find three fasta files in ``FALCON-integrate/FALCON-examples/run/ecoli/data``

2. Create FOFN
~~~~~~~~~~~~~~

Next, create a "file-of-file-names", ("fofn") listing the full path of each fasta file, one per line.

.. code-block:: bash

	/my/path/to/data/ecoli.1.subreads.fasta
	/my/path/to/data/ecoli.2.subreads.fasta
	/my/path/to/data/ecoli.3.subreads.fasta
	

3. Download configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The e-coli config file is available in the :ref:`fc_run.cfg` section.


   
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
	module load python/2.7.9 gcc/4.9.2 git

	# source build
	cd /path/to/my/build/FALCON-integrate
	source env.sh

	# navigate to job directory
	cd /path/to/my/job_dir

	# run it!
	fc_run.py fc_run.cfg


To initiate the FALCON run, I just submit my job to the scheduler with a qsub command:

.. code-block:: bash

	qsub run_falcon.sh



.. _SGE: http://gridscheduler.sourceforge.net/htmlman/manuals.html

Assessing Your Run
------------------

Refer to the pipeline document to understand the FALCON job directory structure, 
sequence of commands, and output files created.

The first step of the FALCON job is to create, then partition of read data into 
