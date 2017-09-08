import falcon_kit.mains.gen_gfa_v1 as mod
import helpers
import pytest
import os

expected_coord_map = {}
expected_coord_map['000000F'] = {
                                    '000092122:B': 0,
                                    '000081654:B': 33726,
                                    '000034462:B': 43849,
                                    '000061403:B': 45201,
                                    '000021348:B': 55125,
                                    '000062240:B': 60959,
                                    '000083779:B': 61821,
                                    '000019819:E': 67383,
                                    '000063672:E': 68767,
                                    '000026565:E': 69240,
                                    '000050047:B': 71411,
                                }
expected_coord_map['000001F'] = {
                                    '000070651:E': 0,
                                    '000018109:E': 10077,
                                    '000068978:E': 13843,
                                    '000010548:B': 18912,
                                    '000100559:E': 16491,
                                    '000006846:B': 21001,
                                    '000065052:B': 39169,
                                    '000071922:E': 41892,
                                    '000076878:E': 44343,
                                    '000000861:E': 45009,
                                    '000001755:B': 60097,
                                }
expected_coord_map['000002F'] = {
                                    '000088930:E': 0,
                                    '000008918:E': 15215,
                                    '000100248:B': 18328,
                                    '000085315:B': 23179,
                                    '000071965:E': 25036,
                                    '000082497:E': 31071,
                                }
expected_coord_map['000003F'] = {
                                    '000084518:E': 0,
                                    '000011674:E': 9432,
                                    '000057445:B': 32528,
                                }
expected_coord_map['000004F'] = {
                                    '000014727:E': 0,
                                    '000024020:E': 5238,
                                    '000060868:E': 8473,
                                }

expected_contig_len = {'000000F': 71411, '000001F': 60097, '000002F': 31071, '000003F': 32528, '000004F': 8473}

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

def test_load_tiling_paths():
    p_ctg_tiling_path_file = os.path.join(helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_path, p_edge_to_ctg = mod.load_tiling_paths(p_ctg_tiling_path_file, 'P')

    assert(sorted(p_path.keys()) == sorted(['000000F', '000001F', '000002F', '000003F', '000004F']))

    for ctg_id, path in p_path.iteritems():
        for edge in path:
            v, w, b, e, l, idt, etype = edge
            assert((v, w) in p_edge_to_ctg)
            assert(p_edge_to_ctg[(v, w)] == (ctg_id, etype))

def test_calc_node_coords():
    # The p_ctg_tiling_path_1 is a normal tiling path file.
    p_ctg_tiling_path_file = os.path.join(helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_paths, p_edge_to_ctg = mod.load_tiling_paths(p_ctg_tiling_path_file, 'P')

    ctg_id = '000000F';
    coord_map, contig_len = mod.calc_node_coords(p_paths[ctg_id])
    shared_items = set(expected_coord_map[ctg_id].items()) & set(coord_map.items())
    assert(len(shared_items) == len(coord_map))
    assert(expected_contig_len[ctg_id] == contig_len)

    # The p_ctg_tiling_path_2 has two degenerative cases:
    # - 000000F which has an inner cycle
    # - 000001F which has an out-of-order edge
    # - 000002F which is circular (this is a valid case)
    p_ctg_tiling_path_file = os.path.join(helpers.get_test_data_dir(), 'p_ctg_tiling_path_2')
    p_paths, p_edge_to_ctg = mod.load_tiling_paths(p_ctg_tiling_path_file, 'P')

    ctg_id = '000000F';
    with pytest.raises(Exception) as e_info:
        coord_map, contig_len = mod.calc_node_coords(p_paths[ctg_id])

    ctg_id = '000001F';
    with pytest.raises(Exception) as e_info:
        coord_map, contig_len = mod.calc_node_coords(p_paths[ctg_id])

    ctg_id = '000002F';
    coord_map, contig_len = mod.calc_node_coords(p_paths[ctg_id])
    assert(contig_len == 18473)

    # Test for an empty tiling path.
    coord_map, contig_len = mod.calc_node_coords([])
    assert(not coord_map)
    assert(contig_len == 0)

def test_calc_tiling_paths_len():
    p_ctg_tiling_path_file = os.path.join(helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_path, p_edge_to_ctg = mod.load_tiling_paths(p_ctg_tiling_path_file, 'P')
    p_coords, p_ctg_len = mod.calc_tiling_paths_len(p_path)

    for ctg_id in p_coords.keys():
        shared_items = set(expected_coord_map[ctg_id].items()) & set(p_coords[ctg_id].items())
        assert(len(shared_items) == len(p_coords[ctg_id]))
        assert(expected_contig_len[ctg_id] == p_ctg_len[ctg_id])

def test_filter_tiling_paths_by_len():
    p_ctg_tiling_path_file = os.path.join(helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_path, p_edge_to_ctg = mod.load_tiling_paths(p_ctg_tiling_path_file, 'P')
    p_coords, p_ctg_len = mod.calc_tiling_paths_len(p_path)

    p_path = mod.filter_tiling_paths_by_len(p_path, p_ctg_len, 0)
    assert(sorted(p_path.keys()) == sorted(['000000F', '000001F', '000002F', '000003F', '000004F']))

    p_path = mod.filter_tiling_paths_by_len(p_path, p_ctg_len, 10000)
    assert(sorted(p_path.keys()) == sorted(['000000F', '000001F', '000002F', '000003F']))

    p_path = mod.filter_tiling_paths_by_len(p_path, p_ctg_len, 35000)
    assert(sorted(p_path.keys()) == sorted(['000000F', '000001F']))

    p_path = mod.filter_tiling_paths_by_len(p_path, p_ctg_len, 100000)
    assert(sorted(p_path.keys()) == sorted([]))

def test_load_seqs():
    fasta_file = os.path.join(helpers.get_test_data_dir(), 't1.fa')

    seqs_to_load = set(['30a5633d_129405_0'])
    load_bases = True
    seq_bases, seq_len = mod.load_seqs(fasta_file, seqs_to_load, load_bases)
    assert(len(seq_bases.keys()) == 1)
    assert(len(seq_len.keys()) == 1)
    for k, seq in seq_bases.iteritems():
        assert(k in seq_len)
        assert(len(seq) == seq_len[k])

    seqs_to_load = set(['30a5633d_129405_0'])
    load_bases = False
    seq_bases, seq_len = mod.load_seqs(fasta_file, seqs_to_load, load_bases)
    assert(len(seq_bases.keys()) == 0)
    assert(len(seq_len.keys()) == 1)
    assert('30a5633d_129405_0' in seq_len)
    assert(seq_len['30a5633d_129405_0'] == 9148)

    seqs_to_load = set([])  # Empty filter set loads everything.
    load_bases = False
    seq_bases, seq_len = mod.load_seqs(fasta_file, seqs_to_load, load_bases)
    assert(len(seq_bases.keys()) == 0)
    assert(len(seq_len.keys()) == 1)
    assert('30a5633d_129405_0' in seq_len)
    assert(seq_len['30a5633d_129405_0'] == 9148)

    seqs_to_load = set([])  # Empty filter set loads everything.
    load_bases = True
    seq_bases, seq_len = mod.load_seqs(fasta_file, seqs_to_load, load_bases)
    assert(len(seq_bases.keys()) == 1)
    assert(len(seq_len.keys()) == 1)
    for k, seq in seq_bases.iteritems():
        assert(k in seq_len)
        assert(len(seq) == seq_len[k])

    seqs_to_load = set(['wrong_header'])
    load_bases = False
    seq_bases, seq_len = mod.load_seqs(fasta_file, seqs_to_load, load_bases)
    assert(len(seq_bases.keys()) == 0)
    assert(len(seq_len.keys()) == 0)

    seqs_to_load = set(['wrong_header'])
    load_bases = True
    seq_bases, seq_len = mod.load_seqs(fasta_file, seqs_to_load, load_bases)
    assert(len(seq_bases.keys()) == 0)
    assert(len(seq_len.keys()) == 0)

    # Use a list instead of set.
    seqs_to_load = ['30a5633d_129405_0']
    load_bases = False
    seq_bases, seq_len = mod.load_seqs(fasta_file, seqs_to_load, load_bases)
    assert(len(seq_bases.keys()) == 0)
    assert(len(seq_len.keys()) == 1)
    assert('30a5633d_129405_0' in seq_len)
    assert(seq_len['30a5633d_129405_0'] == 9148)

def test_format_gfa_path_line():
    pass;

def test_gfa_from_assembly():
    pass;

def test_get_gfa_links_from_tiling_paths():
    pass;

def test_get_gfa_links_from_gfa():
    pass;
