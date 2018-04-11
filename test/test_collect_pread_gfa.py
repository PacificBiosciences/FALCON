import falcon_kit.mains.collect_pread_gfa as mod
import helpers
import pytest
import os

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

def test_load_seqs():
    # load_seqs(fasta_fn, store_only_seq_len)
    pass

def test_load_pread_overlaps():
    # load_pread_overlaps(fp_in)
    pass

def test_load_sg_edges():
    # load_sg_edges(fp_in)
    pass

def test_add_node():
    # add_node(gfa_graph, v, preads_dict)
    pass

def test_add_edge():
    # add_edge(gfa_graph, v, w, edge_split_line, preads_overlap_dict, sg_edges_dict)
    pass

def test_add_tiling_paths_to_gfa():
    # add_tiling_paths_to_gfa(gfa_graph, tiling_paths, preads_dict, preads_overlap_dict, sg_edges_dict)
    pass

def test_add_string_graph_to_gfa():
    # add_string_graph_to_gfa(gfa_graph, sg_edges_list, utg_data, ctg_paths, preads_dict, preads_overlap_dict, sg_edges_dict)
    pass

def test_main_1():
    pass
