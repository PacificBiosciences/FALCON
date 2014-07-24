Quick Note on FALCON Assembly Output Format
============================================

After running `falcon_asm.py`, the following files will be generated

- `edges_list`: the list of edges in the assembled string graph
- `unit_edge_paths`: the path of each unitig
- `unit_edges.dat`: the path and the sequence of each unitig
- `unitgs.fa`: fasta file of all unitigs
- `all_tigs_paths`: the path of all contigs
- `all_tigs.fa`: the sequences of all contigs
- `primary_tigs_paths`: the path of the primary contigs
- `primary_tigs.fa`: the sequences of the initial primary contigs
- `bundle_edges`: the edges and paths of each "string bundles"

After running `falcon_fixasm.py`, it generates the following files

- `primary_tigs_c.fa`: the final primary contigs
- `primary_tigs_paths_c`: the path of the final primary contigs
- `all_tiling_path_c`: the "tiling" path of all contigs
- `primary_tigs_node_pos_c`: the positions of the nodes in each of the primary contigs

The format of each node is the identifier of the DNA fragment followed by `:B` or `:E` indicating the
end of the read that is corresponding to the node.

The `egdes_list` file has a simple 4 column format: `in_node out_node edge_label overlap_length`.
 
Here is an example how edges are represented in the `egdes_list` file:
	
	00099576_1:B 00101043_0:B 00101043_0:1991-0 14333
	00215514_0:E 00025025_0:B 00025025_0:99-0 14948
	00223367_0:E 00146924_0:B 00146924_0:1188-0 8452
	00205542_0:E 00076625_0:B 00076625_0:396-0 11067

The `edge_label`, e.g. `00101043_0:1991-0`, encodes the correspondent sequence of the edge from the DNA fragment. The
edge `00099576_1:B -> 00101043_0:B` has a sequence from read `00101043_0` base 1991 to 0.


The `unit_edge_paths` file contains the information of the path of each unitig. Each line represents 
an unitig. For example, the unitig `00001c` is represented as:

	>00001c-00169881_0:B-00121915_0:E-133 00169881_0:B 00201238_0:E 00137179_0:E 00142410_0:B 
     00223493_0:B 00208425_0:B 00102538_0:E 00160115_0:E  ... 00122905_0:E 00121915_0:E

The full unitig id `00001c-00169881_0:B-00121915_0:E-133` includes the unique serial number `00001c`, the begin node `00169881_0:B` and the end node `00121915_0:E` followed by the number of nodes 133 in the path. The rest of the fields list the full path node by node.

The `primary_tigs_paths` and `all_tigs_paths` have the same format as the `unit_edge_paths` except the edges in the path are the unitig edges rather than the edges in the original string graph.

The `unit_edges.dat` contains not only the begin nodes, the end nodes and the paths of the unitigs but also the full sequences of the unitigs.  It has simple 4 column format `begin node`, `end node`, `path`, `sequence`. The different nodes in the path are delimited by `-`.  

The sequence identifiers in `all_tigs.fa` also encode the relationship between different contigs. For example:

	$ grep ">" all_tigs.fa | head -15
	>0000-0000 2e8a7078_130260_0:B-02eca7b8_135520_0:E
	>0000-0001 6edbcd5c_128868_0:E-3353572d_72448_963:E
	>0000-0002 2f1c350c_15083_0:E-8c92434f_60400_0:E
	>0000-0003 02eca7b8_135520_0:B-02030999_5577_0:B
	>0000-0004-u 53756d78_87035_13099:B-d850f3f2_135807_0:E
	>0000-0005-u 80ae02b0_43730_1168:B-4901e842_5163_2833:B
	>0000-0006-u e1709413_155764_0:E-e55b636f_50757_0:E
	>0000-0007-u e56a70f0_80897_1520:E-06734432_150537_0:E
	>0000-0008-u 1ab64aad_59082_807:E-6f9ad27e_23458_5638:E
	>0000-0009-u 1a88ddf4_21715_0:B-9eb4f7d7_79023_11041:E
	>0000-0010-u ada57c82_24446_0:E-4ce44ebc_41426_0:E
	>0000-0011-u 49704ee2_54679_0:B-a9ced3cc_90191_1410:E
	>0000-0012-u b3728b6f_59022_233:E-bd1579e4_160424_0:B

All these sequences have the same first field `0000`. It means all these contigs are initialized from the same "string bundles". If the second field is `0000`, it means that sequence is the primary contig of this bundle. The rest are the "associated contigs". The second field in the identifier simply indicates the begin and the end node of the contigs.

After running `falcon_fixasm.py`, some of the primary contigs could be broken apart into smaller pieces. For example:
	
	$ grep ">" primary_tigs_c.fa |  head -15
	>0000_00
	>0001_00
	>0001_01
	>0001_02
	>0002_00
	>0002_01

In this case, the initial primary contig `0000` (`0000-0000` in the `all_tigs.fa` file) is intact. However, the `0001-0000` has been broken into 3 primary contigs `0001_00`, `0001_01`, and `0001_02`.

Some of the associated contigs might be caused by sequencing / consensus errors or missing overlapping information. Running `falcon_dedup.py` compares the associated contigs to the corresponding sequences in the primary contigs. If the identity is high, namely not large scale variants found, they will be removed. Mummer3 (Nucmer) package is used and is necessary for this step. `falcon_dedup.py` generates a file called `a_nodup.fa` which contains the non-redundant associate contigs.


Input File Format For FalconSense
---------------------------------

The `falcon_sense.py` generates consensus from a set of raw sequences.

The input is a stream of sequences. Each row has two columns.  Different set of reads are delimited by `+ +` and the file should be ended by `+ +`.  Here is an example

	seed_id1 ACTACATACATACTTA...
	read_id2 TCTGGCAACACTACTTA...
	...
	- -
	seed_id2 ACTACATACATACTTA...
	read_id3 TCTGGCAACACTACTTA...
	...
	- -
	+ +

In this case, if there are enough coverage to correct `seed_id1` and `seed_id2`, the `falcon_sense.py` will generate two consensus sequences (labeled with `seed_id1` and `seed_id2`) in fasta format to `stdout`.

Final Note
----------

1. Typically, the size of `unitgs.fa` will be roughly twice of the genome size, since the file contains both dual edges from each overlap. In the process of the assembly, only one of the dual edges will be used in the final contigs.  

2. The relation between the associate contigs and the primary contigs can be simply identified by the begin and the end nods of the associted contigs. One can easily constructed the corresponding sequences in the primary contigs for identify the variants between them.

3. One can construct a unitig graph from the `unit_edge_paths` files, the graph is typically much smaller than the initial string graph which is more convenient for visualization for understanding the assembly/genome structure.

4. The `-` and `:` characters are used as delimiter for parsing, so the initial reads identifier should not have these two characters. 
