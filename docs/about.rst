.. image:: media/falcon_icon2.png
   :height: 200px
   :width: 200 px
   :alt: FALCON Assembler
   :align: right


.. _about:

About FALCON
============

Overview
--------

``FALCON`` and ``FALCON-Unzip`` are *de novo* genome assemblers for PacBio long reads, also known as 
single-molecule real-time (SMRT) sequences. ``FALCON`` is a diploid-aware assembler 
which follows the hierarchical genome assembly process (HGAP_) and is optimized for 
large genome assembly (e.g. non-microbial). ``FALCON`` produces a set of :term:`primary contigs <primary 
contig>`(p-contigs),
which represent the backbone of the genome sequence, as well as :term:`associate contigs <associated contig>` (a-contigs),
which represent divergent allelic variants. Each a-contig is associated with a homologous
genomic region on an p-contig.

``FALCON-Unzip`` is a true diploid assembler. It takes the contigs from 
``FALCON`` and phases the reads based on heterozygous SNPs identified in the initial 
assembly. It then produces a set of partially phased :term:`primary contigs <primary contig>` and fully phased
:term:`haplotigs <haplotig>` which represent divergent haplotyes.


Detailed Description
--------------------

The hierarchical genome assembly process proceeds in two rounds. The first round of assembly involves the selection of seed reads, 
or the longest reads in the dataset (user-defined :ref:`length_cutoff <length_cutoff>`). All shorter reads are aligned to 
the seed reads, in 
order to generate consensus sequences with high accuracy. We refer to these as pre-assembled reads but they can also be 
thought of as 
“error corrected” reads. During the pre-assembly process, seed reads may be split or trimmed at regions of low read 
coverage (user-defined `min_cov` for :ref:`falcon_sense_option <falcon_sense_option>`). The performance of the pre-assembly 
process is captured in the `pre-assembly stats file.
<http://pb-falcon.readthedocs.io/en/latest/tutorial.html#raw-and-pread-coverage-and-quality>`_

In the next round of HGAP, the :term:`preads <pread>`, are aligned to each other and assembled into 
genomic contigs.

.. image:: media/HGAP.png

.. image:: media/Fig1.png

For more complex genomes assembled with ``FALCON``, 
"bubbles" in the contig-assembly graph that result from structural variation between haplotypes may be resolved as associate 
and primary contigs. The unzip process will extend haplotype phasing beyond "bubble" regions, increasing the amount of phased 
contig sequence. It is important to note that
while individual haplotype blocks are phased, phasing does not extend **between** haplotigs. Thus, in part **C)** of the 
figure above, **haplotig_1** and **haplotig_2** may originate from different parental haplotypes. Additional information is 
needed to phase the haplotype blocks with each other.

Associate contig IDs contain the name of their primary contig but the precise location of alignment must be determined with third party 
tools such as NUCmer_. For example, in a ``FALCON`` assembly, `000123F-010-01` is an associated contig to primary contig 
`000123F`. In a ``FALCON-Unzip`` assembly, `000123F_001` is a haplotig of primary contig `000123F`.

Below are examples of alignments between associate and primary contigs from ``FALCON``, and haplotigs and primary contigs 
from ``FALCON-Unzip``. Alignments were built with NUCmer_ and visualized with Assemblytics_. Precise coordinates 
may be obtained with the show-coords_ utilty from MUMmer_. 

.. image:: media/dotplots.png


Choosing an Assembler: HGAP4 vs FALCON vs FALCON-Unzip 
------------------------------------------------------

HGAP4
~~~~~

We recommend ``HGAP4``, part of the SMRT Link web-based analysis suite, for genomes of known complexity, no larger than 
human (3Gb or 
smaller), 
although underlying 
compute resources for your SMRT Link instance will influence performance and feasibility. The assembly
process for ``HGAP4`` in the SMRT Link GUI (graphical user interface) is identical to ``FALCON`` at the command line, besides 
differences in 
compute resource configuration and minor differences in directory structure. The ``HGAP4`` pipeline by default includes a round of 
genome "polishing" 
which employs the ``resequencing`` pipeline.

``HGAP4`` RESULTS ARE NOT COMPATIBLE WITH ``FALCON-Unzip`` AT THIS TIME!


``HGAP4`` inputs are a PacBio subread BAM_ dataset, either Sequel or RSII. The FASTA and FASTQ files output from ``HGAP4`` are a concatenation of the primary 
and associate contigs, which are output from ``FALCON`` as separate files. 


Command Line
~~~~~~~~~~~~

Users more comfortable at the command line may use ``FALCON`` for genomes of any size 
or complexity. Command line inputs are FASTA files of Sequel or RSII subreads. Command-line ``FALCON`` does not automatically polished the assembly. If a user 
wishes, asembly polishing may 
be run using the ``resequencing`` pipeline of pbsmrtpipe_ (available for command-line installation using the SMRT_Link_ download, see 
SMRT_Tools_Reference_Guide_ for 
installation instructions). Resequencing requires PacBio subread BAM_ inputs.

We recommend the ``FALCON-Unzip`` module for heterozygous or outbred organisms that are diploid or higher ploidy. Users wishing to run 
``FALCON-Unzip`` must do so only after running ``FALCON`` on the 
command line. ``HGAP4`` IS NOT COMPATIBLE WITH ``FALCON-UNZIP``! The ``FALCON-Unzip`` module requires both FASTA and PacBio BAM_ inputs for subreads. 


References
----------

`Chin et al. (2016). Phased diploid genome assembly with single-molecule real-time sequencing. Nature Methods. 13(12), 1050.  
<http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4035.html>`_

`Chin, et al. (2013). Nonhybrid, finished microbial genome assemblies from long-read SMRT sequencing data. Nature Methods. 10(6), 563.
<http://www.nature.com/nmeth/journal/v10/n6/full/nmeth.2474.html>`_


.. _HGAP: http://www.nature.com/nmeth/journal/v10/n6/full/nmeth.2474.html
.. _NUCmer: http://mummer.sourceforge.net/manual/#nucmer
.. _assemblytics: http://qb.cshl.edu/assemblytics/
.. _MUMmer: http://mummer.sourceforge.net/manual/
.. _show-coords: http://mummer.sourceforge.net/manual/#coords
.. _pbsmrtpipe: http://pbsmrtpipe.readthedocs.io/en/master/getting_started.html
.. _SMRT_Link: http://www.pacb.com/support/software-downloads/
.. _SMRT_Tools_Reference_Guide: http://programs.pacificbiosciences.com/l/1652/2017-02-01/3rzxn6/184345/SMRT_Tools_Reference_Guide__v4.0.0_.pdf
.. _BAM: http://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
