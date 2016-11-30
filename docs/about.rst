.. image:: falcon_icon2.png
   :height: 200px
   :width: 200 px
   :alt: Falcon Assembler
   :align: right


.. _about:

About Falcon
============

``FALCON`` and ``FALCON-Unzip`` are *de novo* genome assembly programs for PacBio long reads, also known as 
single-molecule real-time (SMRT) sequences. ``FALCON`` is a diploid-aware assembler 
which follows the hierarchical genome assemble process (HGAP_) and is optimized for 
large genome assembly (e.g. non-microbial). ``FALCON`` produces a set of :term:`primary contigs <primary contig>`,
which represent the backbone of the genome sequence, as well as :term:`associate contigs <associated contig>`
which represent divergent allelic variants. Each a-contig is associated with a homologous
genomic region on an p-contig.

``FALCON-Unzip`` is a true diploid assembler. It takes the genome assembly graph from 
``FALCON`` and phases the reads based on heterozygous SNPs identified in the initial 
assembly. It then produces a set of partially phased :term:`primary contigs <primary contig>` and fully phased
:term:`haplotigs <haplotig>` which represent divergent haplotyes.

.. image:: Fig1.png


References
----------

`Chin et al. Phased diploid genome assembly with single-molecule real-time sequencing. Nat. Meth. 17 Oct 2016
<http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4035.html>`_

.. _HGAP: http://www.nature.com/nmeth/journal/v10/n6/full/nmeth.2474.html



