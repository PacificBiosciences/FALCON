from __future__ import unicode_literals
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

expected_contig_len = {'000000F': 71411, '000001F': 60097,
                       '000002F': 31071, '000003F': 32528, '000004F': 8473}


def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass


def test_main_1(capsys):
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--p-ctg-tiling-path', os.path.join(test_dir, 'p_ctg_tiling_path'),
            '--a-ctg-tiling-path', os.path.join(test_dir, 'a_ctg_tiling_path'),
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--a-ctg-fasta', os.path.join(test_dir, 'a_ctg.fa'),
            '--sg-edges-list', os.path.join(test_dir, 'sg_edges_list'),
            '--utg-data', os.path.join(test_dir, 'utg_data'),
            '--ctg-paths', os.path.join(test_dir, 'ctg_paths'),
            '--add-string-graph',
            '--write-reads',
            '--write-contigs',
            '--min-p-len', '0',
            '--min-a-len', '0',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir,
                                 'expected-1-sg-r-c.gfa')
    helpers.assert_filecmp(out, expected_path)


def test_main_2(capsys):
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--p-ctg-tiling-path', os.path.join(test_dir, 'p_ctg_tiling_path'),
            '--a-ctg-tiling-path', os.path.join(test_dir, 'a_ctg_tiling_path'),
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--a-ctg-fasta', os.path.join(test_dir, 'a_ctg.fa'),
            '--sg-edges-list', os.path.join(test_dir, 'sg_edges_list'),
            '--utg-data', os.path.join(test_dir, 'utg_data'),
            '--ctg-paths', os.path.join(test_dir, 'ctg_paths'),
            # '--add-string-graph',
            '--write-reads',
            '--write-contigs',
            '--min-p-len', '0',
            '--min-a-len', '0',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir,
                                 'expected-2-tiling-r-c.gfa')
    helpers.assert_filecmp(out, expected_path)


def test_main_3(capsys):
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--p-ctg-tiling-path', os.path.join(test_dir, 'p_ctg_tiling_path'),
            '--a-ctg-tiling-path', os.path.join(test_dir, 'a_ctg_tiling_path'),
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--a-ctg-fasta', os.path.join(test_dir, 'a_ctg.fa'),
            '--sg-edges-list', os.path.join(test_dir, 'sg_edges_list'),
            '--utg-data', os.path.join(test_dir, 'utg_data'),
            '--ctg-paths', os.path.join(test_dir, 'ctg_paths'),
            # '--add-string-graph',
            # '--write-reads',
            '--write-contigs',
            '--min-p-len', '0',
            '--min-a-len', '0',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir,
                                 'expected-3-tiling-no_r-c.gfa')
    helpers.assert_filecmp(out, expected_path)


def test_main_4(capsys):
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--p-ctg-tiling-path', os.path.join(test_dir, 'p_ctg_tiling_path'),
            '--a-ctg-tiling-path', os.path.join(test_dir, 'a_ctg_tiling_path'),
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--a-ctg-fasta', os.path.join(test_dir, 'a_ctg.fa'),
            '--sg-edges-list', os.path.join(test_dir, 'sg_edges_list'),
            '--utg-data', os.path.join(test_dir, 'utg_data'),
            '--ctg-paths', os.path.join(test_dir, 'ctg_paths'),
            # '--add-string-graph',
            # '--write-reads',
            # '--write-contigs',
            '--min-p-len', '0',
            '--min-a-len', '0',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir,
                                 'expected-4-tiling-no_r-no_c.gfa')
    helpers.assert_filecmp(out, expected_path)


def test_main_5(capsys):
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--p-ctg-tiling-path', os.path.join(test_dir, 'p_ctg_tiling_path'),
            '--a-ctg-tiling-path', os.path.join(test_dir, 'a_ctg_tiling_path'),
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--a-ctg-fasta', os.path.join(test_dir, 'a_ctg.fa'),
            '--sg-edges-list', os.path.join(test_dir, 'sg_edges_list'),
            '--utg-data', os.path.join(test_dir, 'utg_data'),
            '--ctg-paths', os.path.join(test_dir, 'ctg_paths'),
            '--add-string-graph',
            # '--write-reads',
            # '--write-contigs',
            '--min-p-len', '0',
            '--min-a-len', '0',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir,
                                 'expected-5-sg-no_r-no_c.gfa')
    helpers.assert_filecmp(out, expected_path)


def test_main_6(capsys):
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--p-ctg-tiling-path', os.path.join(test_dir, 'p_ctg_tiling_path'),
            '--a-ctg-tiling-path', os.path.join(test_dir, 'a_ctg_tiling_path'),
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--a-ctg-fasta', os.path.join(test_dir, 'a_ctg.fa'),
            '--sg-edges-list', os.path.join(test_dir, 'sg_edges_list'),
            '--utg-data', os.path.join(test_dir, 'utg_data'),
            '--ctg-paths', os.path.join(test_dir, 'ctg_paths'),
            # '--add-string-graph',
            # '--write-reads',
            # '--write-contigs',
            '--min-p-len', '10000',
            '--min-a-len', '10000',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir,
                                 'expected-6-tiling-no_r-no_c-minlen.gfa')
    helpers.assert_filecmp(out, expected_path)


def test_load_tiling_paths():
    p_ctg_tiling_path_file = os.path.join(
        helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_path, p_edge_to_ctg = mod.load_tiling_paths(p_ctg_tiling_path_file, 'P')

    assert(sorted(p_path.keys()) == sorted(
        ['000000F', '000001F', '000002F', '000003F', '000004F']))

    for ctg_id, path in p_path.iteritems():
        for edge in path:
            v, w, b, e, l, idt, etype = edge
            assert((v, w) in p_edge_to_ctg)
            assert(p_edge_to_ctg[(v, w)] == (ctg_id, etype))


def test_load_tiling_paths_from_stream():
    # This tests a normal case.
    p_ctg_tiling_path_file = os.path.join(
        helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_paths = {}
    edge_to_ctg = {}
    with open(p_ctg_tiling_path_file) as f:
        p_paths, edge_to_ctg = mod.load_tiling_paths_from_stream(f, 'P')
    assert(sorted(p_paths.keys()) == sorted(
        ['000000F', '000001F', '000002F', '000003F', '000004F']))
    for ctg_id, path in p_paths.iteritems():
        for edge in path:
            v, w, b, e, l, idt, etype = edge
            assert((v, w) in edge_to_ctg)
            assert(edge_to_ctg[(v, w)] == (ctg_id, etype))


def test_calc_node_coords():
    # The p_ctg_tiling_path_1 is a normal tiling path file.
    p_ctg_tiling_path_file = os.path.join(
        helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_paths, p_edge_to_ctg = mod.load_tiling_paths(p_ctg_tiling_path_file, 'P')

    ctg_id = '000000F'
    coord_map, contig_len = mod.calc_node_coords(p_paths[ctg_id])
    shared_items = set(expected_coord_map[ctg_id].items()) & set(
        coord_map.items())
    assert(len(shared_items) == len(coord_map))
    assert(expected_contig_len[ctg_id] == contig_len)

    # The p_ctg_tiling_path_2 has two degenerative cases:
    # - 000000F which has an inner cycle
    # - 000001F which has an out-of-order edge
    # - 000002F which is circular (this is a valid case)
    p_ctg_tiling_path_file = os.path.join(
        helpers.get_test_data_dir(), 'p_ctg_tiling_path_2')
    p_paths, p_edge_to_ctg = mod.load_tiling_paths(p_ctg_tiling_path_file, 'P')

    # Allow cycles, but the node's coord gets overwritten.
    ctg_id = '000000F'
    coord_map, contig_len = mod.calc_node_coords(p_paths[ctg_id])
    assert(coord_map['000081654:B'] == 55125)

    # Do not allow unsorted graphs.
    ctg_id = '000001F'
    with pytest.raises(Exception) as e_info:
        coord_map, contig_len = mod.calc_node_coords(p_paths[ctg_id])

    # Allow circular graphs.
    ctg_id = '000002F'
    coord_map, contig_len = mod.calc_node_coords(p_paths[ctg_id])
    assert(contig_len == 18473)

    # Test for an empty tiling path.
    coord_map, contig_len = mod.calc_node_coords([])
    assert(not coord_map)
    assert(contig_len == 0)


def test_calc_tiling_paths_len():
    p_ctg_tiling_path_file = os.path.join(
        helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_path, p_edge_to_ctg = mod.load_tiling_paths(p_ctg_tiling_path_file, 'P')
    p_coords, p_ctg_len = mod.calc_tiling_paths_len(p_path)

    for ctg_id in p_coords.keys():
        shared_items = set(expected_coord_map[ctg_id].items()) & set(
            p_coords[ctg_id].items())
        assert(len(shared_items) == len(p_coords[ctg_id]))
        assert(expected_contig_len[ctg_id] == p_ctg_len[ctg_id])


def test_filter_tiling_paths_by_len():
    p_ctg_tiling_path_file = os.path.join(
        helpers.get_test_data_dir(), 'p_ctg_tiling_path_1')
    p_path, p_edge_to_ctg = mod.load_tiling_paths(p_ctg_tiling_path_file, 'P')
    _, p_ctg_len = mod.calc_tiling_paths_len(p_path)

    p_path_filtered = mod.filter_tiling_paths_by_len(p_path, p_ctg_len, 0)
    assert(sorted(p_path_filtered.keys()) == sorted(
        ['000000F', '000001F', '000002F', '000003F', '000004F']))

    p_path_filtered = mod.filter_tiling_paths_by_len(p_path, p_ctg_len, 10000)
    assert(sorted(p_path_filtered.keys()) == sorted(
        ['000000F', '000001F', '000002F', '000003F']))

    p_path_filtered = mod.filter_tiling_paths_by_len(p_path, p_ctg_len, 35000)
    assert(sorted(p_path_filtered.keys()) == sorted(['000000F', '000001F']))

    p_path_filtered = mod.filter_tiling_paths_by_len(p_path, p_ctg_len, 100000)
    assert(sorted(p_path_filtered.keys()) == sorted([]))

    # Test a degenerate case where there is no length for a particular contig.
    keys = p_ctg_len.keys()
    p_ctg_len_degenerate = {}
    for i in xrange(1, len(keys)):
        p_ctg_len_degenerate[keys[i]] = p_ctg_len[keys[i]]
    with pytest.raises(Exception) as e_info:
        p_path_filtered = mod.filter_tiling_paths_by_len(
            p_path, p_ctg_len_degenerate, 0)


def test_get_filter_tpbc(tmpdir):
    noop = mod.get_filter_tpbc('')
    x = list()
    assert noop(x) is x

    content = """
aaa
bbb
"""
    only_these_contigs = tmpdir.join('only.these')
    only_these_contigs.write(content)
    some = mod.get_filter_tpbc(str(only_these_contigs))
    assert some({'aaa':1, 'bbb-foo':2, 'xxx':3}) == {'aaa':1, 'bbb-foo':2}
