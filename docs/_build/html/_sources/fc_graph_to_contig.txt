.. graph_to_contig:

#####################
fc_graph_to_contig.py
#####################

The final step in the generation of draft contigs is to find a single path for each contig graph and to generate
sequence accordingly. In the case that a contig graph is not a :term:`simple path`, we find the end-to-end path that
has the most overlapped bases. This is called as the :term:`primary contig`. For each :ref:`compound path` within the
graph, if an alternative path different from primary one is possible, we will construct the :ref:`associated contig`.
In the case where the :term:`associated contigs` are induced by sequencing error, the identity of the
alternative contig and the :term:`primary contig` will be high ( > 99% identity most of time).
In the case where there are true structural variations, there are typically bigger differences between the
:term:`associated contigs <associate contig>` and the :term:`primary contigs <primary contig>`.

Essentially, the script ``fc_graph_to_contig`` generates contigs given sequence data and the final assembly graph.
Currently it generates :ref:`primary contigs <primary contig>` as well as all
:ref:`associated contigs <associated contig>` without any filtering. Some post-processing to remove duplicate
:ref:`associated contigs <associated contig>` induced by errors will generally be necessary.
