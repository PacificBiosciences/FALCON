.. _ovlp_to_graph:


###################
fc_ovlp_to_graph.py
###################

Here is the usage information for running fc_ovlp_to_graph.py:

.. code:: bash

    usage: fc_ovlp_to_graph.py [-h] [--min_len MIN_LEN] [--min_idt MIN_IDT]
                               [--lfc]
                               overlap_file

    a example string graph assembler that is desinged for handling diploid genomes

    positional arguments:
      overlap_file       a file that contains the overlap information.

    optional arguments:
      -h, --help         show this help message and exit
      --min_len MIN_LEN  minimum length of the reads to be considered for
                         assembling
      --min_idt MIN_IDT  minimum alignment identity of the reads to be considered
                         for assembling
      --lfc              use local flow constraint method rather than best overlap
                         method to resolve knots in string graph

In some case, you might want to lower the min_idt to keep more overlap or increase min_len to reduce the number of
overlap used for constructing the contig after the overlap filtering step. The ``--lfc`` toggles the rule for resolving
local knots in the graph. If ``--lfc`` is not specified, "the best overlapped edge" will be kept when there are multiple
in- or out- edges from a node while the others will be removed.

The first stage of the assembly is to construct the initial string graph and classify each edges in the string graph.
``sg_edges_list`` contained the information of the information of the edges in the full string graph and the
classification. For example, 5 edges are shown in the five lines of the file below

.. code:: bash

    $ head -5 sg_edges_list
    000017363:B 000007817:E 000007817 10841 28901 10841 99.52 TR
    000015379:E 000004331:B 000004331 6891 0 18178 99.35 TR
    000006813:B 000000681:E 000000681 7609 23795 7616 99.72 TR
    000002258:E 000002505:B 000002505 5850 0 17215 99.62 TR
    000013449:B 000012565:B 000012565 3317 0 20570 99.72 G


The first two columns indicates the in and out node of the edge. The node notation contains two files operated by :.
The first field is the read identifier. The second field is either B or E. B is the ``5'`` end of the read and E is the
``3'`` end of the reads. The next three field indicates the corresponding sequences of the edges. In this example, the
edge in the first line contains the sequence from read ``000007817`` base ``[10841, 28901)``. If the second coordinate is
smaller than the first one, it means the corresponded sequence is reverse complimented. The next two column are the
number of overlapped base and the overlap identity. The final column is the classification. Currently, there are 4
different types ``G``, ``TR``, ``R``, and ``S``. An edge with type ``G`` is used for the final string graph. A ``TR``
means the edge is transitive reducible. ``R`` means the edge is removed during the local repeat resolution and ``S``
means the edge is likely to be a "spur" which only one ends is connected.

The initial string graph is further to be simplified into a set of "unitig" edges. The ``utg_data`` file contains the
details of each unitig. Each line in the file represents a unitig. The first three fields are "start node", "via node",
and "end node". Two untigs might have the same "start node" and "end node", so we need another "via node" to uniquely
identify the unitigs. Here is an example of the utg_data files:


.. code:: bash

    $ head -10 utg_data
    000015696:B 000009028:B 000016941:B contained 16438 134865 000015696:B~000006612:B~000002456:B~000014643:B~000007407:B~000015939:E~000009028:B~000016941:B
    000010623:B 000015633:B 000014991:B contained 30158 18666 000010623:B~000015633:B~000014991:B
    000015636:B 000002245:B 000010757:E contained 15402 40356 000015636:B~000002245:B~000010757:E
    000014184:E NA 000012028:E compound 14895 56928 000014184:E~000012765:E~000012028:E|000014184:E~000007953:B~000012028:E
    000010757:B NA 000015636:E compound 15402 40356 000010757:B~000002245:E~000015636:E|000010757:B~000014783:E~000015636:E
    000014184:E 000007953:B 000012028:E contained 14792 32932 000014184:E~000007953:B~000012028:E
    000010623:B NA 000014991:B compound 30148 163627 000010623:B~000015633:B~000014991:B|000010623:B~000001407:B~000014991:B
    000012028:B 000012765:B 000014184:B contained 19137 56928 000012028:B~000000382:E~000012765:B~000014184:B
    000016941:B 000003353:B 000008783:B simple 88381 615439 000016941:B~000003353:B~000010261:B~000011789:E~000017006:B~000016307:B~...
    000014991:B 000013790:E 000002926:B simple 392373 2274104 000014991:B~000013790:E~000004614:B~000003329:B~000004898:B~000000461:B~000017105:E~...

The forth field indicates the type of the unitigs, the fifth field is the estimate length of the unitig and the six
field is the total number of overlapped bases in the unitig. There are three kinds of unitigs: "simple", "contained",
and "compound". "Simple" unitigs are those unitigs which are just a simple path (every node has one in- and one
out-edge except the begining and ending nodes of the path.) It is represented by a list of nodes which each node
is separated by ~ characters in the 7th column. The "contained" contigs are simple path but those unitigs are also
part of other "compound" paths. The "compound" unitigs represents bubble-like subgraph in the graph. While it is not
"simple", it has well defined in- and out- nodes and they are treated as a single unit when the contigs are
constructed. The structure inside a "compound" unitig can be from biological nature or sequencing/alignment errors.
Each edge in the "compound" unitig subgraph are encoded explicitly as a collection of simple contained unitigs in
the 7th column. The contained unitigs within a compound unitig are separated by the ``|`` character.

The file ctg_paths encodes the graph for each contig after the unitigs are analyzed and put into contigs. Each line
has 7 columns. The first column is the contig ID. The contig ID are just the serial numbers followed by R or F.
Two contigs with same serial number but different endings are "dual" to each other. Namely, they are constructed
from "dual" edges and they are mostly reverse complemented to each other except near the ends of the contigs.
The second column is the type of contig. If a unitig is circular (the beginning node and the ending node are the
same), then it will be marked as ``ctg_circular``. Everything else will be ``ctg_linear``. In some case, even a contig
is marked as ``ctg_linear``, it can be still a circular contig if the beginning node and the ending node are the
same but it is not a "simple" path. One can detect that by checking the beginning and ending nodes if necessary.

The third field indicates the first unitig in the contig in the form of begin_node~via_node~end_node. The fourth
field is the ending node of the contig. The 5th and 6th fields are the estimated length and the overlapped
based of the contig respectively. The final column are the unitigs in the contig. The three node format unitig
IDs are separated by ``|``


