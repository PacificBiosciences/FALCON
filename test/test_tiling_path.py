import falcon_kit.tiling_path as mod
import helpers
import pytest
import os
from StringIO import StringIO

######################################################
# Define an example input test dataset which will be reused.
######################################################
test_1_lines = {
                '000000F': """000000F 000092122:B 000081654:B 000081654 33726 0 16258 99.45
000000F 000081654:B 000034462:B 000034462 10123 0 25619 99.85
000000F 000034462:B 000061403:B 000061403 1352 0 24447 99.96""",
                '000001F': """000001F 000092122:B 000081654:B 000081654 33726 0 16258 99.45""",
                '000002F': """000002F 000014727:E 000024020:E 000024020 15242 20480 15242 99.48
000002F 000024020:E 000060868:E 000060868 19975 23210 19975 99.76"""
}
test_1_expected_coord_map = {
                        '000000F': {'000092122:B': 0, '000081654:B': 33726, '000034462:B': 43849, '000061403:B': 45201},
                        '000001F': {'000092122:B': 0, '000081654:B': 33726},
                        '000002F': {'000014727:E': 0, '000024020:E': 5238, '000060868:E': 8473},
                            }
test_1_expected_contig_len = { '000000F': 45201, '000001F': 33726, '000002F': 8473 }
test_1_input_tp_as_text = '\n'.join([val for key, val in test_1_lines.iteritems()])

######################################################
# Define a new dest dataset for testing placement,
# and also extracting the subpaths.
######################################################
test_placement_1_p_path_as_text = """
000000F 000092122:B 000081654:B 000081654 33726 0 16258 99.45
000000F 000081654:B 000034462:B 000034462 10123 0 25619 99.85
000000F 000034462:B 000061403:B 000061403 1352 0 24447 99.96
000000F 000061403:B 000021348:B 000021348 9924 0 20804 99.74
000000F 000021348:B 000062240:B 000062240 5834 0 27313 99.67
000000F 000062240:B 000083779:B 000083779 862 0 30696 99.79
000000F 000083779:B 000019819:E 000019819 31327 36889 31331 99.94
000000F 000019819:E 000063672:E 000063672 29861 31245 29861 99.96
000000F 000063672:E 000026565:E 000026565 28347 28820 28363 99.87
000000F 000026565:E 000050047:B 000050047 2171 0 23431 99.62
000001F 000092122:B 000081654:B 000081654 33726 0 16258 99.45
000002F 000014727:E 000024020:E 000024020 15242 20480 15242 99.48
000002F 000024020:E 000060868:E 000060868 19975 23210 19975 99.76
"""
test_placement_1_a_path_as_text = """
000000F-1 000034462:B 000070651:E 000070651 0 10000 10000 99.80
000000F-1 000070651:E 000018109:E 000018109 16449 26526 16471 99.68
000000F-1 000018109:E 000068978:E 000068978 24989 28755 24989 99.65
000000F-1 000068978:E 000019819:E 000019819 31327 36889 31331 99.94
"""
test_placement_1_expected = {
    '000000F': {'000000F-1': (43849, 67383, '000000F', '000000F-1', '000034462:B', '000019819:E')},
}
test_placement_1_expected_coords = {
    '000092122:B': 0,       '000081654:B': 33726,   '000034462:B': 43849,
    '000061403:B': 45201,   '000021348:B': 55125,   '000062240:B': 60959,
    '000083779:B': 61821,   '000019819:E': 67383,   '000063672:E': 68767,
    '000026565:E': 69240,   '000050047:B': 71411,
}


def test_tiling_path_1():
    """
    This is a normal test case.
    """

    fp_in = StringIO(test_1_input_tp_as_text)

    # Run test.
    contig_lens = None
    whitelist_seqs = None
    result = mod.load_tiling_paths_from_stream(fp_in, contig_lens=contig_lens, whitelist_seqs=whitelist_seqs)

    # Validate.
    assert(sorted(result.keys()) == ['000000F', '000001F', '000002F'])

    # This checks several methods:
    #   - The entire TilingPathEdge: the values are directly parsed via constructor
    #     called from the TilingPath class. The correctness is tested by callin
    #     the dump_as_split_lines() which should reconstruct the original input line.
    #   - TilingPath has a method `dump_as_split_lines` which returns all such
    #     split lines.
    for key, lines in test_1_lines.iteritems():
        sl = [line.strip().split() for line in lines.splitlines()]
        assert(result[key].dump_as_split_lines() == sl)

    # Check the coordinates.
    assert(result['000000F'].coords == {'000092122:B': 0, '000081654:B': 33726, '000034462:B': 43849, '000061403:B': 45201})
    assert(result['000001F'].coords == {'000092122:B': 0, '000081654:B': 33726})
    assert(result['000002F'].coords == {'000014727:E': 0, '000024020:E': 5238, '000060868:E': 8473})

    assert(result['000000F'].contig_len == 45201)
    assert(result['000001F'].contig_len == 33726)
    assert(result['000002F'].contig_len == 8473)

    assert(result['000000F'].v_to_edge == {'000092122:B': 0, '000081654:B': 1, '000034462:B': 2})
    assert(result['000000F'].w_to_edge == {'000081654:B': 0, '000034462:B': 1, '000061403:B': 2})
    assert(result['000001F'].v_to_edge == {'000092122:B': 0})
    assert(result['000001F'].w_to_edge == {'000081654:B': 0})
    assert(result['000002F'].v_to_edge == {'000014727:E': 0, '000024020:E': 1})
    assert(result['000002F'].w_to_edge == {'000024020:E': 0, '000060868:E': 1})

def test_tiling_path_2():
    """
    This is a normal test case.
    Run with a dict specifying contig lengths. This should offset the coordinates.
    """

    fp_in = StringIO(test_1_input_tp_as_text)

    # Run test.
    contig_lens = {'000000F': 50000, '000001F': 40000, '000002F': 10000}
    whitelist_seqs = None
    result = mod.load_tiling_paths_from_stream(fp_in, contig_lens=contig_lens, whitelist_seqs=whitelist_seqs)

    # Validate.
    assert(sorted(result.keys()) == ['000000F', '000001F', '000002F'])

    for key, lines in test_1_lines.iteritems():
        sl = [line.strip().split() for line in lines.splitlines()]
        assert(result[key].dump_as_split_lines() == sl)

    # Check the coordinates.
    offset = contig_lens['000000F'] - 45201
    assert(result['000000F'].coords == {'000092122:B': 0 + offset, '000081654:B': 33726 + offset, '000034462:B': 43849 + offset, '000061403:B': 45201 + offset})

    offset = contig_lens['000001F'] - 33726
    assert(result['000001F'].coords == {'000092122:B': 0 + offset, '000081654:B': 33726 + offset})

    offset = contig_lens['000002F'] - 8473
    assert(result['000002F'].coords == {'000014727:E': 0 + offset, '000024020:E': 5238 + offset, '000060868:E': 8473 + offset})

    for key, tp in result.iteritems():
        assert(tp.contig_len == contig_lens[key])

    assert(result['000000F'].v_to_edge == {'000092122:B': 0, '000081654:B': 1, '000034462:B': 2})
    assert(result['000000F'].w_to_edge == {'000081654:B': 0, '000034462:B': 1, '000061403:B': 2})
    assert(result['000001F'].v_to_edge == {'000092122:B': 0})
    assert(result['000001F'].w_to_edge == {'000081654:B': 0})
    assert(result['000002F'].v_to_edge == {'000014727:E': 0, '000024020:E': 1})
    assert(result['000002F'].w_to_edge == {'000024020:E': 0, '000060868:E': 1})

def test_tiling_path_3():
    """
    This is a normal test case.
    Test the whitelist filter.
    """

    fp_in = StringIO(test_1_input_tp_as_text)

    # Run test.
    contig_lens = {'000000F': 50000, '000001F': 40000, '000002F': 10000}
    whitelist_seqs = set(['000000F', '000002F', 'key_that_doesnt_exist'])
    result = mod.load_tiling_paths_from_stream(fp_in, contig_lens=contig_lens, whitelist_seqs=whitelist_seqs)

    # Validate.
    assert(sorted(result.keys()) == ['000000F', '000002F'])

    for key in result.keys():
        lines = test_1_lines[key]
        sl = [line.strip().split() for line in lines.splitlines()]
        assert(result[key].dump_as_split_lines() == sl)

    offset = contig_lens['000000F'] - 45201
    assert(result['000000F'].coords == {'000092122:B': 0 + offset, '000081654:B': 33726 + offset, '000034462:B': 43849 + offset, '000061403:B': 45201 + offset})

    offset = contig_lens['000002F'] - 8473
    assert(result['000002F'].coords == {'000014727:E': 0 + offset, '000024020:E': 5238 + offset, '000060868:E': 8473 + offset})

    for key, tp in result.iteritems():
        assert(tp.contig_len == contig_lens[key])

    assert(result['000000F'].v_to_edge == {'000092122:B': 0, '000081654:B': 1, '000034462:B': 2})
    assert(result['000000F'].w_to_edge == {'000081654:B': 0, '000034462:B': 1, '000061403:B': 2})
    assert(result['000002F'].v_to_edge == {'000014727:E': 0, '000024020:E': 1})
    assert(result['000002F'].w_to_edge == {'000024020:E': 0, '000060868:E': 1})

def test_tiling_path_4_expect_crash():
    """
    This test attempts to create a tiling path from edges which are out of order.
    An exception should be raised when initializing the TilingPath object.
    The `calc_node_coords` should raise an exception for node '000021348:B'.
    """

    test_3_lines = {}
    test_3_lines['000000F'] = test_1_lines['000000F'] + """\n000000F 000062240:B 000083779:B 000083779 862 0 30696 99.79"""
    test_3_lines['000001F'] = test_1_lines['000001F']
    test_3_lines['000002F'] = test_1_lines['000002F']
    test_3_input_tp_as_text = '\n'.join([val for key, val in test_3_lines.iteritems()])

    fp_in = StringIO(test_3_input_tp_as_text)

    contig_lens = None
    whitelist_seqs = None

    with pytest.raises(Exception):
        result = mod.load_tiling_paths_from_stream(fp_in, contig_lens=contig_lens, whitelist_seqs=whitelist_seqs)

def test_tiling_path_5(tmpdir):
    """
    This is a normal test case for loading from file. The results should be the
    same as loading from stream.
    """

    contig_lens = None
    whitelist_seqs = None

    tp_file = tmpdir.join('test_tiling_path')
    tp_file.write(test_1_input_tp_as_text)
    result = mod.load_tiling_paths(str(tp_file), contig_lens=contig_lens, whitelist_seqs=whitelist_seqs)

    fp_in = StringIO(test_1_input_tp_as_text)
    expected = mod.load_tiling_paths_from_stream(fp_in, contig_lens=contig_lens, whitelist_seqs=whitelist_seqs)

    assert(result.keys() == expected.keys())

    for key in result:
        result_tp = result[key]
        expected_tp = expected[key]
        assert(result_tp.coords == expected_tp.coords)
        assert(result_tp.dump_as_split_lines() == expected_tp.dump_as_split_lines())
        assert(result_tp.contig_len == expected_tp.contig_len)
        assert(result_tp.v_to_edge == expected_tp.v_to_edge)
        assert(result_tp.w_to_edge == expected_tp.w_to_edge)

def test_calc_node_coords_1():
    """
    Run a normal test.
    """

    fp_in = StringIO(test_1_input_tp_as_text)

    # Run test.
    # The calc_node_coords requires a list of TilingEdge objects.
    # Parsing of tiling paths and loading of TilingEdge objects was already
    # tested above, se here we just reuse the mechanism to get tp.edges.

    expected_coord_map = {}
    expected_coord_map['000000F'] = {'000092122:B': 0, '000081654:B': 33726, '000034462:B': 43849, '000061403:B': 45201}
    expected_coord_map['000001F'] = {'000092122:B': 0, '000081654:B': 33726}
    expected_coord_map['000002F'] = {'000014727:E': 0, '000024020:E': 5238, '000060868:E': 8473}

    expected_contig_len = {}
    expected_contig_len['000000F'] = 45201
    expected_contig_len['000001F'] = 33726
    expected_contig_len['000002F'] = 8473

    contig_lens = expected_contig_len
    whitelist_seqs = None
    tps = mod.load_tiling_paths_from_stream(fp_in, contig_lens=contig_lens, whitelist_seqs=whitelist_seqs)

    for key, tp in tps.iteritems():
        coord_map, contig_len = mod.calc_node_coords(tp.edges, first_node_offset=0)
        assert(coord_map == expected_coord_map[key])
        assert(contig_len == expected_contig_len[key])

def test_calc_node_coords_2():
    """
    Test on empty input.
    """
    coord_map, contig_len = mod.calc_node_coords([], first_node_offset=0)

    assert(len(coord_map.keys()) == 0)
    assert(contig_len == 0)

def test_calc_node_coords_3():
    """
    Run a normal test.
    """

    fp_in = StringIO(test_1_input_tp_as_text)

    # Run test.
    # The calc_node_coords requires a list of TilingEdge objects.
    # Parsing of tiling paths and loading of TilingEdge objects was already
    # tested above, se here we just reuse the mechanism to get tp.edges.

    first_node_offset = {}
    first_node_offset['000000F'] = 50000 - 45201
    first_node_offset['000001F'] = 40000 - 33726
    first_node_offset['000002F'] = 10000 - 8473

    expected_contig_len = {}
    expected_contig_len['000000F'] = 50000 # 45201
    expected_contig_len['000001F'] = 40000 # 33726
    expected_contig_len['000002F'] = 10000 # 8473

    expected_coord_map = {}
    offset = first_node_offset['000000F']
    expected_coord_map['000000F'] = {'000092122:B': 0 + offset, '000081654:B': 33726 + offset, '000034462:B': 43849 + offset, '000061403:B': 45201 + offset}
    offset = first_node_offset['000001F']
    expected_coord_map['000001F'] = {'000092122:B': 0 + offset, '000081654:B': 33726 + offset}
    offset = first_node_offset['000002F']
    expected_coord_map['000002F'] = {'000014727:E': 0 + offset, '000024020:E': 5238 + offset, '000060868:E': 8473 + offset}

    contig_lens = expected_contig_len
    whitelist_seqs = None
    tps = mod.load_tiling_paths_from_stream(fp_in, contig_lens=contig_lens, whitelist_seqs=whitelist_seqs)

    for key, tp in tps.iteritems():
        coord_map, contig_len = mod.calc_node_coords(tp.edges, first_node_offset=first_node_offset[key])
        assert(coord_map == expected_coord_map[key])
        assert(contig_len == expected_contig_len[key])

def test_calc_node_coords_1():
    """
    Attempt to run on edges which are out of order.
    This should raise an exception.
    """

    test_input_tp_as_text = test_1_lines['000000F'] + """\n000000F 000062240:B 000083779:B 000083779 862 0 30696 99.79"""
    all_sl = [line.strip().split() for line in test_input_tp_as_text.splitlines()]
    edges = [mod.TilingPathEdge(sl) for sl in all_sl]
    with pytest.raises(Exception):
        coord_map, contig_len = mod.calc_node_coords(edges, first_node_offset=0)

def test_find_a_ctg_placement_1():
    """
    Normal test case.
    The find_a_ctg_placement method expects a valid primary contig tiling path dict and a
    valid associate contig tiling path dict as inputs.
    Validation of the construction of TilingPath objects is performed in the tests
    above; here we take these objects for granted.
    """
    fp_in = StringIO(test_placement_1_p_path_as_text)
    p_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)
    fp_in = StringIO(test_placement_1_a_path_as_text)
    a_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)

    result = mod.find_a_ctg_placement(p_paths, a_paths)

    for key, placement in result.iteritems():
        assert(placement == test_placement_1_expected[key])

def test_find_a_ctg_placement_2():
    """
    Test empty a_ctg paths.
    """
    fp_in = StringIO(test_placement_1_p_path_as_text)
    p_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)
    a_paths = {}

    result = mod.find_a_ctg_placement(p_paths, a_paths)

    assert(len(result.keys()) == 0)

def test_find_a_ctg_placement_3():
    """
    Test empty p_ctg paths. If the primary contig cannot be found, this should
    throw an exception.
    """
    p_paths = {}
    fp_in = StringIO(test_placement_1_a_path_as_text)
    a_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)

    with pytest.raises(Exception):
        result = mod.find_a_ctg_placement(p_paths, a_paths)

def test_get_subpath_1():
    """
    Check the case when the beginning and the end of the contig are provided.
    This should extract the entire path.
    """

    fp_in = StringIO(test_placement_1_p_path_as_text)
    p_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)

    p_path = p_paths['000000F']

    # The result is composed of (subpath, start_in_subpath, end_in_subpath).
    test_start = 0
    test_end = p_path.contig_len
    result = p_path.get_subpath(test_start, test_end)

    expected = [val.get_split_line() for val in p_path.edges], test_start, test_end

    assert(result == expected)

def test_get_subpath_2():
    """
    Test with a start coord < 0, which is possible if the contig was improper.
    """

    fp_in = StringIO(test_placement_1_p_path_as_text)
    p_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)

    p_path = p_paths['000000F']

    # The result is composed of (subpath, start_in_subpath, end_in_subpath).
    test_start = -100
    test_end = p_path.contig_len
    result = p_path.get_subpath(test_start, test_end)
    subpath, start_in_subpath, end_in_subpath = result

    expected = [val.get_split_line() for val in p_path.edges], test_start, test_end

    assert(result == expected)

def test_get_subpath_3():
    """
    Test the end coordinate larger than the contig length. This is also possible
    if the contig was improper, but is frowned upon.
    """

    fp_in = StringIO(test_placement_1_p_path_as_text)
    p_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)

    p_path = p_paths['000000F']

    # The result is composed of (subpath, start_in_subpath, end_in_subpath).
    test_start = -100
    test_end = p_path.contig_len + 1000
    result = p_path.get_subpath(test_start, test_end)
    subpath, start_in_subpath, end_in_subpath = result

    expected = [val.get_split_line() for val in p_path.edges], test_start, test_end

    assert(result == expected)

def test_get_subpath_4():
    """
    Normal test - extract a subpath between internal coordinates.
    """

    fp_in = StringIO(test_placement_1_p_path_as_text)
    p_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)
    p_path = p_paths['000000F']

    # The result is composed of (subpath, start_in_subpath, end_in_subpath).
    test_start = 43849 + 10 # Start is within the third edge.
    test_end = 68767 + 100  # End is within the 9th edge
    result = p_path.get_subpath(test_start, test_end)

    subpath, start_in_subpath, end_in_subpath = result

    expected = [val.get_split_line() for val in p_path.edges[2:9]], 10, (test_end - 43849)

    assert(result == expected)

def test_get_subpath_5():
    """
    Test edge case, when the selected coordinates are right at the ends of
    the tiling path edges.
    """

    fp_in = StringIO(test_placement_1_p_path_as_text)
    p_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)
    p_path = p_paths['000000F']

    # The result is composed of (subpath, start_in_subpath, end_in_subpath).
    test_start = 43849 + 0 # Start is within the third edge.
    test_end = 69240 + 0  # End is within the 9th edge
    result = p_path.get_subpath(test_start, test_end)

    subpath, start_in_subpath, end_in_subpath = result

    expected = [val.get_split_line() for val in p_path.edges[2:9]], 0, (test_end - 43849)

    assert(result == expected)

def test_get_subpath_6():
    """
    The end coordinate is only one base into a new edge. This edge should entirely
    be output.
    """

    fp_in = StringIO(test_placement_1_p_path_as_text)
    p_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)
    p_path = p_paths['000000F']

    # The result is composed of (subpath, start_in_subpath, end_in_subpath).
    test_start = 43849 + 0 # Start is within the third edge.
    test_end = 69240 + 1  # End is within the 9th edge
    result = p_path.get_subpath(test_start, test_end)

    subpath, start_in_subpath, end_in_subpath = result

    expected = [val.get_split_line() for val in p_path.edges[2:10]], 0, (test_end - 43849)

    assert(result == expected)

def test_get_subpath_7():
    """
    When both coordinates are <= 0, then only the first edge should be output.
    """

    fp_in = StringIO(test_placement_1_p_path_as_text)
    p_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)
    p_path = p_paths['000000F']

    # The result is composed of (subpath, start_in_subpath, end_in_subpath).
    test_start = -100
    test_end = -10
    result = p_path.get_subpath(test_start, test_end)

    subpath, start_in_subpath, end_in_subpath = result

    expected = [val.get_split_line() for val in p_path.edges[0:1]], test_start, test_end

    assert(result == expected)

def test_get_subpath_8():
    """
    The end coordinate should not be <= start coordinate.
    """

    fp_in = StringIO(test_placement_1_p_path_as_text)
    p_paths = mod.load_tiling_paths_from_stream(fp_in, contig_lens=None, whitelist_seqs=None)
    p_path = p_paths['000000F']

    test_start = 10
    test_end = 7
    with pytest.raises(Exception):
        result = p_path.get_subpath(test_start, test_end)
