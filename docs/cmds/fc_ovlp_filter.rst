.. ovlp_filter:

#################
fc_ovlp_filter.py
#################


.. code:: bash

    usage: fc_ovlp_filter [-h] [--n_core N_CORE] [--fofn FOFN] [--db DB_FN]
                          [--max_diff MAX_DIFF] [--max_cov MAX_COV]
                          [--min_cov MIN_COV] [--min_len MIN_LEN] [--bestn BESTN]
                          [--stream] [--debug] [--silent]

    a simple multi-processes LAS ovelap data filter

    optional arguments:
      -h, --help           show this help message and exit
      --n_core N_CORE      number of processes used for generating consensus; 0
                           for main process only (default: 4)
      --fofn FOFN          file contains the path of all LAS file to be processed
                           in parallel (default: None)
      --db DB_FN           read db file path (default: None)
      --max_diff MAX_DIFF  max difference of 5' and 3' coverage (default: None)
      --max_cov MAX_COV    max coverage of 5' or 3' coverage (default: None)
      --min_cov MIN_COV    min coverage of 5' or 3' coverage (default: None)
      --min_len MIN_LEN    min length of the reads (default: 2500)
      --bestn BESTN        output at least best n overlaps on 5' or 3' ends if
                           possible (default: 10)
      --stream             stream from LA4Falcon, instead of slurping all at once;
                           can save memory for large data (default: False)
      --debug, -g          single-threaded, plus other aids to debugging (default:
                           False)
      --silent             suppress cmd reporting on stderr (default: False)


Not all overlaps are "independent", so it is possible to impose some filtering step to reduce computation and
assembly graph complexity. For example, if a read is fully contained in another read, the overlap information
between these two reads does not provide extra information for re-constructing the genome. Also, due to the
transitive property of the overlapping relations, a lot of overlap information can be simply inferred. In fact,
the first stage for constructing contigs are to remove the transitive reducible edges. It means that we might
just needs the "best n overlaps" in the ``5'`` or ``3'`` ends. The ``--bestn`` parameter in :ref:`overlap_filtering_setting`
option can be used to control the maximum overlap reported for each read.

Another useful heuristics is to only keep reads that have average ``5'`` and ``3'`` coverage. That's because if a read
ends in a repeat, it might have higher than normal coverage at the end which is a repeat. And such reads do not
provide much value for uniquely resolving the related repeats. We can filter them out and hopefully there are
reads which span through the repeats and have "normal" unique anchors on both ends. Also, if the coverage is too
low on one end of a read, it could be just too many errors or sequencing artifacts over there. Such reads create
"spurs" in the assembly graph which are typically filtered out anyway. The --max_cov and --min_cov are used for
filtering reads that have too high or too low overlaps.

The filtering scripts also allows filtering out some "split" reads. If a read have very unequal coverage between
the ``5'`` and ``3'`` ends, it can be also a signal that one end is a repeat. The ``--max_diff`` parameter can be used to
filter out the reads where one ends has much more coverage than the other end.

What is the right numbers used for these parameters? These parameters may the most tricky ones to be set right.
If the overall coverage of the error corrected reads longer than the length cut off is known and reasonable high
(e.g. greater than 20x), it might be safe to set ``--min_cov`` to be 5, max_cov to be three times of the average
coverage and the max_diff to be twice of the average coverage. However, in low coverage case, it might better
to set ``--min_cov`` to be one or two. A helper script called :doc:`fc_ovlp_stats`` can help to dump the number of the
``3'`` and ``5'`` overlap of a given length cutoff, you can plot the distribution of the number of overlaps to make a
better decision.

One can also set the ``--max_diff`` and ``--max_cov`` to be really high to avoid any filtering if that is preferred
in some cases.

This filtering process will certainly filter out information about high copy repeats. Namely, those repeats will
likely to be filtered out totally and do not appear in the final assembly. If you are interested in those repeats
even though they may not be able to placed within some longer contig, you will probably want to avoid filtering
them out or process them differently. In general, it might be more efficient and useful to process those repeats
separately. Including them in the assembly process typically does not help much for getting better contiguity and
maybe messy for post-processing with current algorithms. I think it is a very interesting but also very challenging
bioinformatics topic on how to process these repeats better for improving assembly beside understand the nature
of these repeats.