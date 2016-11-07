.. _commands:

Commands
========

Dazzler commands
----------------

These commands are part of Gene Meyer's Dazzler Suite of tools `Dazzler Blog <http://dazzlerblog.wordpress.com>`_
FALCON relies on a slightly modified version of Gene Meyer's code that can be found
`here <https://github.com/cschin/DALIGNER>`_


.. _daligner:

`daligner <https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide>`_:
    ``daligner`` is controlled by :ref:`pa_HPCdaligner_option` and :ref:`ovlp_HPCdaligner_option`.

    To limit memory, one can use the ``-M`` option. For human assembly, we've tested with ``-M 32`` for using 32G RAM for
    each daligner. Other possibilities are under investigation.

    For more details on daligner options, see the `Dazzler Blog <http://dazzlerblog.wordpress.com>`


.. _DB2Fasta:

`DB2fasta <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    info

.. _DBdump:

`DBdump <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    info

.. _DBdust:

`DBdust <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    stuff

.. _DBsplit:

`DBsplit <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    The total number of jobs that are run is determined by how one "splits" the sequence database. You should read
    Gene Myers's blog `Dazzler Blog <http://dazzlerblog.wordpress.com>` carefully to understand how the tuning options,
    :ref:`pa_DBsplit_option` and :ref:`pa_HPCdaligner_option` work. Generally, for large genomes, you should use
    ``-s400`` (400Mb sequence per block) in :ref:`pa_DBsplit_option`. This will make a smaller number of jobs but each
    job will run longer. However, if you have a job scheduler which limits how long a job can run, it might be
    desirable to have a smaller number for the ``-s`` option.

.. _DBstats:

`DBstats <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    info

.. _fasta2DB:
`fasta2DB <https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide>`_:
    Convert a fasta to a dazzler DB. See dazzDB_ command guide for details

.. _HPC.daligner:

`HPC.daligner <https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide>`_:
    info

.. _LA4Falcon:

`LA4Falcon <cmds/LA4Falcon>`_:
    Output data from a Dazzler DB into fasta format for FALCON. You can supply the argument ``-H`` with an integer value
    to filter reads below a given threshold.

.. _LAcheck:

`LAcheck <https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide>`_:
    Check integrity of alignment files.

.. _LAmerge:

`LAmerge <https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide>`_:
    The total number of jobs that are run is determined by how one "splits" the sequence database. You should read
    Gene Myers's blog ( http://dazzlerblog.wordpress.com ) carefully to know how to tune the option pa_DBsplit_option
    and pa_HPCdaligner_option. Generally, for large genomes, you should use -s400 (400Mb sequence per block) in
    pa_DBsplit_option. This will make a smaller number of jobs but each job will run longer. However, if you have a job
    queue system which limits how long a job can run, it might be desirable to have a smaller number for the -s option.

.. _LAsort:

`LAsort <https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide>`_:
    Sort alignment files


FALCON Commands
---------------

.. _DB2Falcon:

:doc:`DB2Falcon <cmds/DB2Falcon>`
    Used to dump dazzler preads.db into FASTA format for subsequent :term:`String Graph` assembly

.. _fc_run:

:doc:`fc_run <cmds/fc_run>`
    This script drives the entire assembly process

.. _fc_consensus:

:doc:`fc_consensus <cmds/fc_consensus>`
    ``fc_consensus`` has many options. You can use the parameter :ref:`falcon_sense_option` to control it.
    In most cases, the ``--min_cov`` and ``--max_n_read`` are the most important options. ``--min_cov`` controls
    when a seed read gets trimmed or broken due to low coverage. ``--max_n_read`` puts a cap on the number of reads
    used for error correction. In highly repetitive genome, you will need to make the value for ``--max_n_read``
    smaller to make sure the consensus code does not waste time aligning repeats. The longest proper overlaps are used
    for correction to reduce the probability of collapsed repeats.

.. _fc_dedup_a_tigs:

:doc:`fc_dedup_a_tigs <cmds/fc_dedup_a_tigs>`
    info

.. _fc_graph_to_contig:

:doc:`fc_graph_to_contig <cmds/fc_graph_to_contig>`
    Generate contigs based on assembly graph

.. _fc_ovlp_to_graph:

:doc:`fc_ovlp_to_graph <cmds/fc_ovlp_to_graph>`
    Generate an assembly graph given a list of overlapping preads.

.. _fc_ovlp_filter:

:doc:`fc_ovlp_filter <cmds/fc_ovlp_to_graph>`
    Filter overlaps based on given criteria

