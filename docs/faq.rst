.. _faq:

Frequently Asked Questions
==========================



Can I start from corrected reads?
+++++++++++++++++++++++++++++++++

Yes. The option :ref:`input_type` can be set to either ``raw`` or ``preads``. In the case of the latter, ``fc_run.py``
will assume the fasta files in :ref:`input_fofn` are all error-corrected reads and it will ignore any error correction
step and go directly into the final assembly overlapping step.

How do I select a length cutoff?
++++++++++++++++++++++++++++++++

The option :ref:`length_cutoff` controls the read length cutoff used during the error correction process and
:ref:`length_cutoff_pr` controls the cutoff used for the final assembly overlapping steps. In the final assembly,
more reads may not lead to a better assembly due to the fact that some of the reads can be noisy and create false
links in the assembly graph. Sometimes you might want to re-run the final steps of the assembly pipeline in
``2-asm-falcon`` with different values for ``--min_len`` in ``run_falcon_asm.sub.sh`` as this step is quick relative
to the overlap detection steps in the earlier stages of the pipeline.

If you're not sure, and you are not compute resource limited, one strategy is to choose a smaller :ref:`length_cutoff`
and do the computation once. Later, we can use a different :ref:`length_cutoff_pr` for getting better assembly.

In general we recommend that you tune the cutoff so that you're left with roughly 15x to 20x for final genome assembly.
If you set :ref:`length_cutoff` equal to ``-1``, FALCON will attempt to autocalculate this cutoff for you.