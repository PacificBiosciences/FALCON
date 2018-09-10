.. _commands:

.. caution:: These documents refer to an obsolete way of installing and running FALCON. They will remain up for historical context and for individuals still using the older version of FALCON/FALCON_unzip.

.. attention:: The current PacBio Assembly suite documentation which includes new bioconda instructions for installing FALCON, FALCON_unzip and their associated dependencies can be found here `pb_assembly <http://github.com/PacificBiosciences/pb-assembly>`_


Commands
========

FALCON Commands
---------------

.. _DB2dalcon:

:doc:`DB2Falcon <cmds/DB2Falcon>`
    Used to dump dazzler preads.db into FASTA format for subsequent :term:`String Graph` assembly

.. _fc_run.py:

:doc:`fc_run.py <cmds/fc_run>`
    This script drives the entire assembly process

.. _fc_consensus.py:

:doc:`fc_consensus.py <cmds/fc_consensus>`
    ``fc_consensus`` has many options. You can use the parameter :ref:`falcon_sense_option` to control it.
    In most cases, the ``--min_cov`` and ``--max_n_read`` are the most important options. ``--min_cov`` controls
    when a seed read gets trimmed or broken due to low coverage. ``--max_n_read`` puts a cap on the number of reads
    used for error correction. In highly repetitive genome, you will need to make the value for ``--max_n_read``
    smaller to make sure the consensus code does not waste time aligning repeats. The longest proper overlaps are used
    for correction to reduce the probability of collapsed repeats.

.. _fc_dedup_a_tigs.py:

:doc:`fc_dedup_a_tigs.py <cmds/fc_dedup_a_tigs>`
    remove duplicated :term:`associated contigs<associated contig>`, mostly induced by tandem repeat alignment
    uncertainty

.. _fc_graph_to_contig.py:

:doc:`fc_graph_to_contig.py <cmds/fc_graph_to_contig>`
    Generate contigs based on assembly graph

.. _fc_ovlp_to_graph.py:

:doc:`fc_ovlp_to_graph.py <cmds/fc_ovlp_to_graph>`
    Generate an assembly graph given a list of overlapping preads.

.. _fc_ovlp_filter.py:

:doc:`fc_ovlp_filter.py <cmds/fc_ovlp_to_graph>`
    Filter overlaps based on given criteria

FALCON_unzip commands
---------------------

:doc:`fc_get_read_hctg_map.py <cmds/fc_get_read_hctg_map>`
    Generate a read-to-contig map

.. _fc_dedup_h_tigs.py:

:doc:`fc_dedup_h_tigs.py <cmds/fc_dedup_h_tigs>`

.. _fc_graphs_to_h_tigs.py:

:doc:`fc_graphs_to_h_tigs.py <cmds/fc_graphs_to_h_tigs>`

.. _fc_ovlp_filter_with_phase.py:

:doc:`fc_ovlp_filter_with_phase.py <cmds/fc_ovlp_filter_with_phase>`

.. _fc_phased_ovlp_to_graph.py:

:doc:`fc_phased_ovlp_to_graph.py <cmds/fc_phased_ovlp_to_graph>`

.. _fc_phasing.py:

:doc:`fc_phasing.py <cmds/fc_phasing>`

.. _fc_phasing_readmap.py:

:doc:`fc_phasing_readmap.py <cmds/fc_phasing_readmap>`

.. _fc_quiver.py:

:doc:`fc_quiver.py <cmds/fc_quiver>`

.. _fc_rr_hctg_track.py:

:doc:`fc_rr_hctg_track.py <cmds/fc_rr_hctg_track>`

.. _fc_select_reads_from_bam.py:

:doc:`fc_select_reads_from_bam.py <cmds/fc_select_reads_from_bam>`

.. _fc_track_reads_htigs.py:

:doc:`fc_track_reads_htigs.py <cmds/fc_track_reads_htigs0>`

.. _fc_unzip.py:

:doc:`fc_unzip.py <cmds/fc_unzip>`

Dazzler commands
----------------

These commands are part of Gene Meyer's Dazzler Suite of tools `Dazzler Blog <http://dazzlerblog.wordpress.com>`_

FALCON relies on a slightly modified version of Gene Meyer's code that can be found
`here <https://github.com/cschin/DALIGNER>`_, but is also bundled with the
`FALCON-integrate <https://github.com/PacificBiosciences/FALCON-integrate.git>`_ github repository.

.. _dazzdaligner:

`daligner <https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide>`_:
    Compare subject sequences to target sequences
    ``daligner`` is controlled by :ref:`pa_HPCdaligner_option <pa_HPCdaligner_option>` and
    :ref:`ovlp_HPCdaligner_option <ovlp_HPCdaligner_option>`.

    To limit memory, one can use the ``-M`` option. For human assembly, we've tested with ``-M 32`` for using 32G RAM for
    each daligner. Other possibilities are under investigation.

    For more details on daligner options, see the `Dazzler Blog <http://dazzlerblog.wordpress.com>`_

.. _dazzDB2fasta:

`DB2fasta <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    The set of .fasta files for the given DB are recreated from the DB exactly as they were input.

.. _dazzDBdump:

`DBdump <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    Like DBshow, DBdump allows one to display a subset of the reads in the DB and select which information to show
    about them including any mask tracks.

.. _dazzDBdust:

`DBdust <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    Runs the symmetric DUST algorithm over the reads in the untrimmed DB

.. _dazzDBsplit:

`DBsplit <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    The total number of jobs that are run is determined by how one "splits" the sequence database. You should read
    Gene Myers's blog `Dazzler Blog <http://dazzlerblog.wordpress.com>` carefully to understand how the tuning options,
    :ref:`pa_DBsplit_option <pa_DBsplit_option>` and :ref:`pa_HPCdaligner_option <pa_HPCdaligner_option>` work. Generally, for large genomes, you should use
    ``-s400`` (400Mb sequence per block) in :ref:`pa_DBsplit_option <pa_DBsplit_option>`. This will make a smaller number of jobs but each
    job will run longer. However, if you have a job scheduler which limits how long a job can run, it might be
    desirable to have a smaller number for the ``-s`` option.

.. _dazzDBstats:

`DBstats <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    Show overview statistics for all the reads in the trimmed data base <path>.db

.. _dazzfasta2DB:

`fasta2DB <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    Convert a fasta to a dazzler DB.

.. _dazzHPC.daligner:

`HPC.daligner <https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide>`_:
    Generates overlap script to run all necessary daligner, LAsort and LAmerge commands

.. _dazzLA4Falcon:

`LA4Falcon <cmds/LA4Falcon>`_:
    Output data from a Dazzler DB into fasta format for FALCON. You can supply the argument ``-H`` with an integer value
    to filter reads below a given threshold.

.. _dazzLAcheck:

`LAcheck <https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide>`_:
    Check integrity of alignment files.

.. _dazzLAmerge:

`LAmerge <https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide>`_:
    Merge the .las files <parts> into a singled sorted file

.. _dazzLAsort:

`LAsort <https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide>`_:
    Sort alignment files
