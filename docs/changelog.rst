.. _changelog:

Changelog
=========


.. _3122018tarball:

3/12/2018
---------


.. _542018tarball:

5/4/2018
--------

FALCON Updates:
+++++++++++++++

Repeat Masking

- Integration of DAMASKER: Tandem repeat masking (done) and general repeat masking (in progress)

Improved default settings for microbial assembly

- Use one longest read per ZMW: reduced chimerism, coverage bias
- Retuned parameters to increase contiguity

New! GFA-2 support

- Assembly graphs now written in both GFA-1 and GFA-2 formats
- Placement coordinates of associate contigs now available in a new "contig.gfa2" file

Performance Improvements

- General workflow and resource specification improvements
- Easier integration of future features with Pbsmrtpipe


FALCON-Unzip Updates:
+++++++++++++++++++++

Improved Haplotig Extraction

- Algorithm and data structure improvements reduce haplotype switching and improve extraction (resolved nested and overlapping haplotigs)
- Can now handle circular contigs!

New! Placement Files

- Haplotig placement (PAF format) generated after Unzip
- Easier integration with FALCON-Phase

Performance Improvements

- Use of Minimap2 instead of BLASR for phasing in Unzip reduces time requirements
- Significantly reduced memory consumption of the final stage of Unzip (preads no longer have to be loaded in memory)
- Unzipping and polishing are now combined in the same workflow and run consecutively.

