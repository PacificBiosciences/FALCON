.. _tutorial:

.. image:: falcon_icon2.png
   :height: 200px
   :width: 200 px
   :alt: Falcon Assembler
   :align: right

Tutorial
========

In this section we will work through a concrete example of running Falcon on real data. 
We will work through the commands and results as well as give you ideas of how to assess 
the perfomance of falcon on your dataset so you can modify parameters and trouble-shoot effectively.

Input Files
-----------

We will use the *E. coli* example dataset provided with the falcon-integrate distribution. Follow the :ref:`quick_start` 
to install Falcon and set up your environment. Next, navigate to the FALCON-integrate/FALCON-examples directory in 
install the dataset:

.. code-block:: bash

	cd FALCON-example
	make run-ecoli
	
This prints to screen the raw data in fasta format as well as some additional commands. 
I'll leave it to you to format the data into pure fasta files. The data can be in 1 or many
fasta files.

Next, create a "file-of-file-names", ("fofn") listing the full path of each fasta file, one per line.

You can use the e-coli config file 

   
Running Falcon
------------

I send all of my Falcon jobs to the scheduler for ease of tracking




Introduction
____________


