.. _fc_consensus:


###############
fc_consensus.py
###############

.. code-block:: bash

    usage: fc_consensus.py [-h] [--n_core N_CORE] [--min_cov MIN_COV]
                           [--min_cov_aln MIN_COV_ALN] [--max_cov_aln MAX_COV_ALN]
                           [--min_len_aln MIN_LEN_ALN] [--min_n_read MIN_N_READ]
                           [--max_n_read MAX_N_READ] [--trim] [--output_full]
                           [--output_multi] [--min_idt MIN_IDT]
                           [--edge_tolerance EDGE_TOLERANCE]
                           [--trim_size TRIM_SIZE]

    a simple multi-processor consensus sequence generator

    optional arguments:
      -h, --help            show this help message and exit
      --n_core N_CORE       number of processes used for generating consensus; 0
                            for main process only (default: 24)
      --min_cov MIN_COV     minimum coverage to break the consensus (default: 6)
      --min_cov_aln MIN_COV_ALN
                            minimum coverage of alignment data; a seed read with
                            less than MIN_COV_ALN average depth of coverage will
                            be completely ignored (default: 10)
      --max_cov_aln MAX_COV_ALN
                            maximum coverage of alignment data; a seed read with
                            more than MAX_COV_ALN average depth of coverage of the
                            longest alignments will be capped, excess shorter
                            alignments will be ignored (default: 0)
      --min_len_aln MIN_LEN_ALN
                            minimum length of a sequence in an alignment to be
                            used in consensus; any shorter sequence will be
                            completely ignored (default: 0)
      --min_n_read MIN_N_READ
                            1 + minimum number of reads used in generating the
                            consensus; a seed read with fewer alignments will be
                            completely ignored (default: 10)
      --max_n_read MAX_N_READ
                            1 + maximum number of reads used in generating the
                            consensus (default: 500)
      --trim                trim the input sequence with k-mer spare dynamic
                            programming to find the mapped range (default: False)
      --output_full         output uncorrected regions too (default: False)
      --output_multi        output multi correct regions (default: False)
      --min_idt MIN_IDT     minimum identity of the alignments used for correction
                            (default: 0.7)
      --edge_tolerance EDGE_TOLERANCE
                            for trimming, the there is unaligned edge leng >
                            edge_tolerance, ignore the read (default: 1000)
      --trim_size TRIM_SIZE
                            the size for triming both ends from initial sparse
                            aligned region (default: 50)
