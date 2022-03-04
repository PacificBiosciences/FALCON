.. _fc_get_read_hctg_map

#######################
fc_get_read_hctg_map.py
#######################

.. code-block:: bash

    usage: fc_get_read_hctg_map.py [-h] [--basedir BASEDIR]

    generate `3-unzip/read_maps/read_to_contig_map` that contains the information
    from the chain of mapping: (contig id, last col) -> (internal p-read id) ->
    (internal raw-read id) -> (original read id) it assumes the 2-asm-
    falcon/read_maps/raw_read_ids and pread_ids are already generated

    optional arguments:
      -h, --help         show this help message and exit
      --basedir BASEDIR  the base working dir of a FALCON assembly (default: ./)