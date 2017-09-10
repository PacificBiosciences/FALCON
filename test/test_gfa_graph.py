import os
import pytest
import helpers
from StringIO import StringIO
import networkx as nx
import falcon_kit.gfa_graph as mod
import falcon_kit.mains.gen_gfa_v1 as gen_gfa_v1
from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.mains.ovlp_to_graph import reverse_end

def test_gfa_graph():
    gfa_graph = mod.GFAGraph()

def test_add_tiling_path():
    # Load the tiling path. These methods are tested in test_gen_gfa_v1.py.
    p_ctg_tiling_path_file = os.path.join(helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_paths, p_edge_to_ctg = gen_gfa_v1.load_tiling_paths(p_ctg_tiling_path_file, 'P')

    # Create a new GFA graph.
    gfa_graph = mod.GFAGraph()

    # Add the tiling paths.
    for ctg_id, path in p_paths.iteritems():
        gfa_graph.add_tiling_path(path, ctg_id)

    # Check if we have the correct number of tiling paths.
    assert(len(gfa_graph.paths.keys()) == len(p_paths.keys()))

    # They should be same as loaded.
    for ctg_id, path in p_paths.iteritems():
        assert(ctg_id in gfa_graph.paths)
        assert(gfa_graph.paths[ctg_id] == path)

def test_add_asm_graph():
    # Load the assembly graph.
    sg_edges_list = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'sg_edges_list')
    utg_data = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'utg_data')
    ctg_paths = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'ctg_paths')
    asm_graph = AsmGraph(sg_edges_list, utg_data, ctg_paths)

    # Add the graph to GFA.
    gfa_graph = mod.GFAGraph()
    gfa_graph.add_asm_graph(asm_graph)

    assert(len(gfa_graph.paths.keys()) == 0)

    expected = {
                ('000000016:B', '000000027:B'): ['000000016:B', '000000027:B', '*', 1540, 99.94, 449, 0, None, None, None, None],
                ('000000005:B', '000000016:B'): ['000000005:B', '000000016:B', '*', 1487, 99.93, 502, 0, None, None, None, None],
                ('000000016:B', '000000025:B'): ['000000016:B', '000000025:B', '*', 1540, 99.94, 449, 0, None, None, None, None],
                ('000000007:B', '000000005:B'): ['000000007:B', '000000005:B', '*', 1980, 99.95, 9, 0, None, None, None, None],
                ('000000018:B', '000000004:B'): ['000000018:B', '000000004:B', '*', 1963, 99.95, 26, 0, None, None, None, None],
                ('000000025:B', '000000018:B'): ['000000025:B', '000000018:B', '*', 1978, 99.95, 11, 0, None, None, None, None]
               }

    assert(len(gfa_graph.edges.keys()) == len(expected.keys()))

    for key, edge in gfa_graph.edges.iteritems():
        assert(key in expected)
        assert(expected[key] == edge)

def test_add_nx_string_graph():
    # Load the assembly graph.
    sg_edges_list = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'sg_edges_list')
    utg_data = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'utg_data')
    ctg_paths = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'ctg_paths')
    asm_graph = AsmGraph(sg_edges_list, utg_data, ctg_paths)

    # The following block is taken from Unzip, graphs_to_h_tigs.py.
    nx_sg = nx.DiGraph()
    arid_to_phase = {}
    for ctg_id in asm_graph.ctg_data.keys():
        ctg_G = asm_graph.get_sg_for_ctg(ctg_id)
        ctg_nodes = set(ctg_G.nodes())
        for v, w in ctg_G.edges():
            vrid = v[:9]
            wrid = w[:9]
            edge_data = asm_graph.sg_edges[ (v, w) ]
            if edge_data[-1] != "G":
                continue

            vphase = arid_to_phase.get(vrid, (-1,0))
            wphase = arid_to_phase.get(wrid, (-1,0))
            if vphase[0] == wphase[0] and vphase[1] != wphase[1]:
                cross_phase = "Y"
            else:
                cross_phase = "N"

            nx_sg.add_node( v, label= "%d_%d" % vphase,
                            phase="%d_%d" % vphase,
                            src="P" )

            nx_sg.add_node( w, label= "%d_%d" % wphase,
                            phase="%d_%d" % wphase,
                            src="P" )

            nx_sg.add_edge(v, w, src="OP", cross_phase = cross_phase)

            # we need to add the complimentary edges as the ctg_graph does not contain the dual edges
            rv = reverse_end(v)
            rw = reverse_end(w)
            nx_sg.add_node( rv, label= "%d_%d" % vphase,
                            phase="%d_%d" % vphase,
                            src="P" )
            nx_sg.add_node( rw, label= "%d_%d" % wphase,
                            phase="%d_%d" % wphase,
                            src="P" )
            nx_sg.add_edge(rw, rv, src="OP", cross_phase = cross_phase)

    # Add the string graph to the GFA.
    gfa_graph = mod.GFAGraph()
    gexf_file = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'sg.gexf')
    nx_sg = nx.read_gexf(gexf_file)
    gfa_graph.add_nx_string_graph(nx_sg)
    # nx.write_gexf(nx_sg, gexf_file)

def wrap_write_gfa_v1_test(use_sg, use_nx, use_tp, write_reads, write_contigs, min_p_len, min_a_len, expected_path):
    # Create a GFA graph.
    gfa_graph = mod.GFAGraph()

    if use_sg:
        # Load the assembly graph.
        sg_edges_list = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'sg_edges_list')
        utg_data = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'utg_data')
        ctg_paths = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'ctg_paths')
        asm_graph = AsmGraph(sg_edges_list, utg_data, ctg_paths)
        # Add the string graph to the GFA.
        gfa_graph.add_asm_graph(asm_graph)

    if use_tp:
        # Load the p_ctg tiling paths.
        p_ctg_tiling_path_file = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'p_ctg_tiling_path')
        p_paths, p_edge_to_ctg = gen_gfa_v1.load_tiling_paths(p_ctg_tiling_path_file, 'P')
        # Add the tiling paths to the GFA.
        for ctg_id, path in p_paths.iteritems():
            _, contig_len = gen_gfa_v1.calc_node_coords(path)
            if contig_len >= min_p_len:
                gfa_graph.add_tiling_path(path, ctg_id)
        a_ctg_tiling_path_file = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'a_ctg_tiling_path')
        a_paths, a_edge_to_ctg = gen_gfa_v1.load_tiling_paths(a_ctg_tiling_path_file, 'P')
        # Add the tiling paths to the GFA.
        for ctg_id, path in a_paths.iteritems():
            _, contig_len = gen_gfa_v1.calc_node_coords(path)
            if contig_len >= min_a_len:
                gfa_graph.add_tiling_path(path, ctg_id)

    if use_nx:
        gexf_file = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'sg.gexf')
        nx_sg = nx.read_gexf(gexf_file)
        gfa_graph.add_nx_string_graph(nx_sg)

    # Init paths to other input files.
    preads_file = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'preads4falcon.fasta')
    p_ctg_fasta = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'p_ctg.fa')
    a_ctg_fasta = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'a_ctg.fa')

    fp_out = StringIO()
    # Run the unit under test.
    gfa_graph.write_gfa_v1(fp_out, preads_file, [p_ctg_fasta, a_ctg_fasta], write_reads, write_contigs)

    # Compare results.
    result = fp_out.getvalue()
    result = result.splitlines()
    expected = [line.strip() for line in open(expected_path).readlines()]
    assert(result == expected)

def test_write_gfa_v1_1():
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    # Test various combinations of options.
    wrap_write_gfa_v1_test(True, False, True, True, True, 0, 0, os.path.join(test_dir, 'expected-1-sg-r-c.gfa'))
    wrap_write_gfa_v1_test(False, False, True, True, True, 0, 0, os.path.join(test_dir, 'expected-2-tiling-r-c.gfa'))
    wrap_write_gfa_v1_test(False, False, True, False, True, 0, 0, os.path.join(test_dir, 'expected-3-tiling-no_r-c.gfa'))
    wrap_write_gfa_v1_test(False, False, True, False, False, 0, 0, os.path.join(test_dir, 'expected-4-tiling-no_r-no_c.gfa'))
    wrap_write_gfa_v1_test(True, False, True, False, False, 0, 0, os.path.join(test_dir, 'expected-5-sg-no_r-no_c.gfa'))
    wrap_write_gfa_v1_test(False, False, True, False, False, 10000, 10000, os.path.join(test_dir, 'expected-6-tiling-no_r-no_c-minlen.gfa'))
    wrap_write_gfa_v1_test(False, True, False, False, False, 0, 0, os.path.join(test_dir, 'expected-7-nx-no_r-no_c.gfa'))
    wrap_write_gfa_v1_test(False, True, True, False, False, 0, 0, os.path.join(test_dir, 'expected-8-nx-tiling-no_r-no_c.gfa'))
    wrap_write_gfa_v1_test(False, True, True, True, True, 0, 0, os.path.join(test_dir, 'expected-9-nx-tiling-r-c.gfa'))

def test_write_gfa_v1_2():
    # Tests a case where a node is added to the graph, but
    # there is no corresponding pread in preads4falcon.fasta file.

    # Create a GFA graph.
    gfa_graph = mod.GFAGraph()

    # Load the p_ctg tiling paths.
    p_ctg_tiling_path_file = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'p_ctg_tiling_path')
    p_paths, p_edge_to_ctg = gen_gfa_v1.load_tiling_paths(p_ctg_tiling_path_file, 'P')
    # Add the tiling paths to the GFA.
    for ctg_id, path in p_paths.iteritems():
        gfa_graph.add_tiling_path(path, ctg_id)

    # Init paths to other input files.
    preads_file = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'preads4falcon.fasta')
    p_ctg_fasta = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'p_ctg.fa')
    a_ctg_fasta = os.path.join(helpers.get_test_data_dir(), 'gfa-1', 'a_ctg.fa')

    write_reads = False
    write_contigs = False

    fp_out = StringIO()

    # Add a node which does not exist in the preads4falcon.fasta file.
    gfa_graph.add_read_from_node('12345:B')

    # Run the unit under test.
    with pytest.raises(Exception) as e_info:
        gfa_graph.write_gfa_v1(fp_out, preads_file, [p_ctg_fasta, a_ctg_fasta], write_reads, write_contigs)

def test_add_read_from_node():
    gfa_graph = mod.GFAGraph()

    gfa_graph.add_read_from_node('123:B')
    assert(len(gfa_graph.read_in_graph) == 1)
    assert('123' in gfa_graph.read_in_graph)

    gfa_graph.add_read_from_node('456:')
    assert(len(gfa_graph.read_in_graph) == 2)
    assert('123' in gfa_graph.read_in_graph)
    assert('456' in gfa_graph.read_in_graph)

    gfa_graph.add_read_from_node('123:B')
    assert(len(gfa_graph.read_in_graph) == 2)
    assert('123' in gfa_graph.read_in_graph)
    assert('456' in gfa_graph.read_in_graph)

    gfa_graph.add_read_from_node('123:E')
    assert(len(gfa_graph.read_in_graph) == 2)
    assert('123' in gfa_graph.read_in_graph)
    assert('456' in gfa_graph.read_in_graph)

    with pytest.raises(Exception) as e_info:
        gfa_graph.add_read_from_node('123')
    with pytest.raises(Exception) as e_info:
        gfa_graph.add_read_from_node(None)
    with pytest.raises(Exception) as e_info:
        gfa_graph.add_read_from_node('')

def test_add_edge():
    gfa_graph = mod.GFAGraph()

    v, w, cigar = '123:B', '456:E', '*'
    overlap_len, overlap_idt, overlap_begin, overlap_end = 10000, 99.9, 0, 9000
    cross_phase, src_graph, ctg_name, type_ = 'N', 'OP', '000000F', 'P'
    gfa_graph.add_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)
    assert(len(gfa_graph.read_in_graph) == 2)
    assert(len(gfa_graph.edges.keys()) == 1)
    assert((v, w) in gfa_graph.edges)

    # Check that multiedges cannot be added.
    gfa_graph.add_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)
    assert(len(gfa_graph.read_in_graph) == 2)
    assert(len(gfa_graph.edges.keys()) == 1)
    assert((v, w) in gfa_graph.edges)

    assert(v.split(':')[0] in gfa_graph.read_in_graph)
    assert(w.split(':')[0] in gfa_graph.read_in_graph)

def test_update_edge():
    gfa_graph = mod.GFAGraph()

    # First, add an edge.
    v, w, cigar = '123:B', '456:E', '*'
    overlap_len, overlap_idt, overlap_begin, overlap_end = 10000, 99.9, 0, 9000
    cross_phase, src_graph, ctg_name, type_ = None, None, None, None
    gfa_graph.add_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)
    assert(len(gfa_graph.edges.keys()) == 1)
    assert((v, w) in gfa_graph.edges)
    assert(gfa_graph.edges[(v, w)] == ['123:B', '456:E', '*', 10000, 99.9, 0, 9000, None, None, None, None])

    # Update the None values and check if they changed.
    v, w, cigar = '123:B', '456:E', '*'
    overlap_len, overlap_idt, overlap_begin, overlap_end = 10000, 99.9, 0, 9000
    cross_phase, src_graph, ctg_name, type_ = 'N', 'OP', '000000F', 'P'
    gfa_graph.update_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)
    assert(len(gfa_graph.edges.keys()) == 1)
    assert((v, w) in gfa_graph.edges)
    assert(gfa_graph.edges[(v, w)] == ['123:B', '456:E', '*', 10000, 99.9, 0, 9000, 'N', 'OP', '000000F', 'P'])

    # Add a new edge
    v, w, cigar = '456:B', '789:E', '*'
    overlap_len, overlap_idt, overlap_begin, overlap_end = 10000, 99.9, 0, 9000
    cross_phase, src_graph, ctg_name, type_ = None, None, None, None
    gfa_graph.add_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)
    assert(len(gfa_graph.edges.keys()) == 2)
    assert((v, w) in gfa_graph.edges)
    assert(gfa_graph.edges[(v, w)] == ['456:B', '789:E', '*', 10000, 99.9, 0, 9000, None, None, None, None])

    # Update the values, but check that non-None values remained the same.
    v, w, cigar = '456:B', '789:E', '*'
    overlap_len, overlap_idt, overlap_begin, overlap_end = 10001, 99.8, 0, 9001
    cross_phase, src_graph, ctg_name, type_ = 'N', 'OP', '000000F', 'P'
    gfa_graph.update_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)
    assert(len(gfa_graph.edges.keys()) == 2)
    assert((v, w) in gfa_graph.edges)
    assert(gfa_graph.edges[(v, w)] == ['456:B', '789:E', '*', 10000, 99.9, 0, 9000, 'N', 'OP', '000000F', 'P'])

    # Degenerate case, update an edge which does not exist.
    v, w, cigar = '4567:B', '789:E', '*'
    overlap_len, overlap_idt, overlap_begin, overlap_end = 10001, 99.8, 0, 9001
    cross_phase, src_graph, ctg_name, type_ = 'N', 'OP', '000000F', 'P'
    with pytest.raises(Exception) as e_info:
        gfa_graph.update_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)

def test_add_or_update_edge():
    gfa_graph = mod.GFAGraph()

    # First, add an edge.
    v, w, cigar = '123:B', '456:E', '*'
    overlap_len, overlap_idt, overlap_begin, overlap_end = 10000, 99.9, 0, 9000
    cross_phase, src_graph, ctg_name, type_ = None, None, None, None
    gfa_graph.add_or_update_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)
    assert(len(gfa_graph.edges.keys()) == 1)
    assert((v, w) in gfa_graph.edges)
    assert(gfa_graph.edges[(v, w)] == ['123:B', '456:E', '*', 10000, 99.9, 0, 9000, None, None, None, None])

    # Update the None values and check if they changed.
    v, w, cigar = '123:B', '456:E', '*'
    overlap_len, overlap_idt, overlap_begin, overlap_end = 10000, 99.9, 0, 9000
    cross_phase, src_graph, ctg_name, type_ = 'N', 'OP', '000000F', 'P'
    gfa_graph.add_or_update_edge(v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_)
    assert(len(gfa_graph.edges.keys()) == 1)
    assert((v, w) in gfa_graph.edges)
    assert(gfa_graph.edges[(v, w)] == ['123:B', '456:E', '*', 10000, 99.9, 0, 9000, 'N', 'OP', '000000F', 'P'])

def test_format_gfa_v1_link_line():
    gfa_graph = mod.GFAGraph()

    # Test an edge with None information.
    v, w, cigar = '123:B', '456:E', '*'
    overlap_len, overlap_idt, overlap_begin, overlap_end = 10000, 99.9, 0, 9000
    cross_phase, src_graph, ctg_name, type_ = None, None, None, None
    edge = [v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_]
    result = gfa_graph.format_gfa_v1_link_line(edge)
    expected = 'L\t123\t-\t456\t+\t*\tol:i:10000\toi:f:99.9\toe:i:9000\tci:Z:NA-NA'
    assert(result == expected)

    # Test an edge with full information.
    v, w, cigar = '456:B', '789:E', '*'
    overlap_len, overlap_idt, overlap_begin, overlap_end = 10000, 99.9, 0, 9000
    cross_phase, src_graph, ctg_name, type_ = 'N', 'OP', '000000F', 'P'
    edge = [v, w, cigar, overlap_len, overlap_idt, overlap_begin, overlap_end, cross_phase, src_graph, ctg_name, type_]
    result = gfa_graph.format_gfa_v1_link_line(edge)
    expected = 'L\t456\t-\t789\t+\t*\tol:i:10000\toi:f:99.9\toe:i:9000\tsg:Z:OP\tcp:Z:N\tci:Z:000000F-P'
    assert(result == expected)

def test_format_gfa_v1_path_line():
    gfa_graph = mod.GFAGraph()

    # Load tiling paths from file.
    p_ctg_tiling_path_file = os.path.join(helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_paths, p_edge_to_ctg = gen_gfa_v1.load_tiling_paths(p_ctg_tiling_path_file, 'P')

    # If seq_len_map is None, all CIGAR operations should be '*'.
    expected = {
                '000000F': 'P\t000000F\t000092122-,000081654-,000034462-,000061403-,000021348-,000062240-,000083779-,000019819+,000063672+,000026565+,000050047-\t*,*,*,*,*,*,*,*,*,*,*',
                '000001F': 'P\t000001F\t000070651+,000018109+,000068978+,000100559+,000010548-,000006846-,000065052-,000071922+,000076878+,000000861+,000001755-\t*,*,*,*,*,*,*,*,*,*,*',
                '000002F': 'P\t000002F\t000088930+,000008918+,000100248-,000085315-,000071965+,000082497+\t*,*,*,*,*,*',
                '000003F': 'P\t000003F\t000084518+,000011674+,000057445-\t*,*,*',
                '000004F': 'P\t000004F\t000014727+,000024020+,000060868+\t*,*,*',
               }
    seq_len_map = None
    for ctg_id, path in p_paths.iteritems():
        path_line = gfa_graph.format_gfa_v1_path_line(ctg_id, path, seq_len_map)
        assert(path_line == expected[ctg_id])

    # The seq_len_map dict is only used for the first read in the path,
    # because it needs to be included completely. The other CIGAR operations
    # are determined directly from the edges.
    expected = {
                '000000F': 'P\t000000F\t000092122-,000081654-,000034462-,000061403-,000021348-,000062240-,000083779-,000019819+,000063672+,000026565+,000050047-\t10000M,33726M,10123M,1352M,9924M,5834M,862M,5562M,1384M,473M,2171M',
                '000001F': 'P\t000001F\t000070651+,000018109+,000068978+,000100559+,000010548-,000006846-,000065052-,000071922+,000076878+,000000861+,000001755-\t10000M,10077M,3766M,2648M,2421M,2089M,18168M,2723M,2451M,666M,15088M',
                '000002F': 'P\t000002F\t000088930+,000008918+,000100248-,000085315-,000071965+,000082497+\t10000M,15215M,3113M,4851M,1857M,6035M',
                '000003F': 'P\t000003F\t000084518+,000011674+,000057445-\t10000M,9432M,23096M',
                '000004F': 'P\t000004F\t000014727+,000024020+,000060868+\t10000M,5238M,3235M',
               }
    for ctg_id, path in p_paths.iteritems():
        # Initialize all reads to a fixed value, just to be safe.
        seq_len_map = {}
        for edge in path:
            v, w = edge[0], edge[1]
            seq_len_map[v.split(':')[0]] = 10000
            seq_len_map[w.split(':')[0]] = 10000
        path_line = gfa_graph.format_gfa_v1_path_line(ctg_id, path, seq_len_map)
        assert(path_line == expected[ctg_id])

    # Test a degenerate case where path is None.
    path_line = gfa_graph.format_gfa_v1_path_line('', None, None)
    assert(path_line == '')
