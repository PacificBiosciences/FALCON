from StringIO import StringIO
import os

import pytest
import networkx as nx

from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.mains.ovlp_to_graph import reverse_end
from falcon_kit.util import system
import falcon_kit.gfa_graph as mod
import falcon_kit.mains.gen_gfa_v1 as gen_gfa_v1
import helpers


def test_gfa_graph():
    gfa_graph = mod.GFAGraph()

def test_add_node_1():
    """
    Test normal usage.
    """
    gfa_graph = mod.GFAGraph()

    gfa_graph.add_node('node1', 4, 'ACTG', tags={}, labels={})
    gfa_graph.add_node('node2', 1000, '*', tags={}, labels={})

    assert(len(gfa_graph.nodes) == 2)

def test_add_node_2():
    """
    Tests that exceptions get raised if parameters are not correct.
    """

    gfa_graph = mod.GFAGraph()

    with pytest.raises(Exception):
        gfa_graph.add_node('', 4, 'ACTG', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_node('node1', -1, 'ACTG', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_node('node1', 4, '', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_node('node1', 4, 'ACTG', tags=[], labels={})
    with pytest.raises(Exception):
        gfa_graph.add_node('node1', 4, 'ACTG', tags={}, labels=[])

def test_add_edge_1():
    """
    Test normal usage.
    """

    gfa_graph = mod.GFAGraph()

    gfa_graph.add_node('node1', 4, 'ACTG', tags={}, labels={})
    gfa_graph.add_node('node2', 1000, '*', tags={}, labels={})

    edge_name = 'edge1'
    source, source_orient = 'node1', '+'
    sink, sink_orient = 'node2', '+'
    source_start, source_end = 4, 4
    sink_start, sink_end = 0, 0
    cigar = '*'

    gfa_graph.add_edge(edge_name, source, source_orient, sink, sink_orient, source_start, source_end, sink_start, sink_end, cigar, tags={}, labels={})

    assert(len(gfa_graph.edges.keys()) == 1)

def test_add_edge_2():
    """
    Tests that exceptions get raised if parameters are not correct.
    """

    gfa_graph = mod.GFAGraph()

    with pytest.raises(Exception):
        gfa_graph.add_edge('', 'node1', '+', 'node2', '+', 4, 4, 0, 0, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', '', '+', 'node2', '+', 4, 4, 0, 0, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', '', '+', 4, 4, 0, 0, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '1', 'node2', '+', 4, 4, 0, 0, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', 'node2', 'z', 4, 4, 0, 0, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', -1, 4, 0, 0, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, -1, 0, 0, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 4, -1, 0, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 4, 0, -1, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 4, 0, 0, '', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 3, 0, 0, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 4, 5, 0, '*', tags={}, labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 4, 0, 0, '*', tags=[], labels={})
    with pytest.raises(Exception):
        gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 4, 0, 0, '*', tags={}, labels=[])

def test_add_path_1():
    """
    Test normal usage.
    """

    gfa_graph = mod.GFAGraph()

    path_nodes = ['node1', 'node2', 'node3', 'node4', 'node5']
    path_cigars = ['500M', '400M', '300M', '200M', '100M']

    gfa_graph.add_path('000000F', path_nodes, path_cigars)

    assert(len(gfa_graph.paths.keys()) == 1)

def test_add_path_2():
    """
    Tests that exceptions get raised if parameters are not correct.
    """

    gfa_graph = mod.GFAGraph()

    with pytest.raises(Exception):
        path_nodes = ['node1', 'node2', 'node3', 'node4', 'node5']
        path_cigars = ['500M', '400M', '300M', '200M', '100M']
        gfa_graph.add_path('', path_nodes, path_cigars)
    with pytest.raises(Exception):
        path_nodes = ['node1']
        path_cigars = ['500M', '400M', '300M', '200M', '100M']
        gfa_graph.add_path('000000F', path_nodes, path_cigars)
    with pytest.raises(Exception):
        path_nodes = ['node1', 'node2', 'node3', 'node4', 'node5']
        path_cigars = ['500M']
        gfa_graph.add_path('000000F', path_nodes, path_cigars)
    with pytest.raises(Exception):
        path_nodes = []
        path_cigars = ['500M']
        gfa_graph.add_path('000000F', path_nodes, path_cigars)
    with pytest.raises(Exception):
        path_nodes = ['node1']
        path_cigars = []
        gfa_graph.add_path('000000F', path_nodes, path_cigars)
    with pytest.raises(Exception):
        path_nodes = ['node1', 'node2', 'node3', 'node4', 'node5']
        path_cigars = ['500M', '400M', '300M', '200M', '100M']
        gfa_graph.add_path('000000F', path_nodes, path_cigars, tags=[], labels={})
    with pytest.raises(Exception):
        path_nodes = ['node1', 'node2', 'node3', 'node4', 'node5']
        path_cigars = ['500M', '400M', '300M', '200M', '100M']
        gfa_graph.add_path('000000F', path_nodes, path_cigars, tags={}, labels=[])

def test_write_gfa_v1():
    gfa_graph = mod.GFAGraph()
    gfa_graph.add_node('node1', 7, 'ACTGAAA', tags={}, labels={})
    gfa_graph.add_node('node2', 10, 'AAACCCGGGT', tags={}, labels={})
    gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 7, 0, 3, '*', tags={}, labels={})
    gfa_graph.add_path('000000F', ['node1', 'node2'], ['4M', '7M'], tags={}, labels={})

    fp_out = StringIO()
    gfa_graph.write_gfa_v1(fp_out)

    result = fp_out.getvalue()
    expected = """H	VN:Z:1.0
S	node1	ACTGAAA	LN:i:7
S	node2	AAACCCGGGT	LN:i:10
L	node1	+	node2	+	3M
P	000000F	node1,node2	4M,7M
"""
    assert(result == expected)

def test_write_gfa_v2():
    gfa_graph = mod.GFAGraph()
    gfa_graph.add_node('node1', 7, 'ACTGAAA', tags={}, labels={})
    gfa_graph.add_node('node2', 10, 'AAACCCGGGT', tags={}, labels={})
    gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 7, 0, 3, '*', tags={}, labels={})
    gfa_graph.add_path('000000F', ['node1', 'node2'], ['4M', '7M'], tags={}, labels={})

    fp_out = StringIO()
    gfa_graph.write_gfa_v2(fp_out)

    result = fp_out.getvalue()
    expected = """H	VN:Z:2.0
S	node1	7	ACTGAAA
S	node2	10	AAACCCGGGT
E	edge1	node1+	node2+	4	7$	0	3	*
"""

    assert(result == expected)

def test_serialize():
    gfa_graph = mod.GFAGraph()
    gfa_graph.add_node('node1', 7, 'ACTGAAA', tags={}, labels={})
    gfa_graph.add_node('node2', 10, 'AAACCCGGGT', tags={}, labels={})
    gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 7, 0, 3, '*', tags={}, labels={})
    gfa_graph.add_path('000000F', ['node1', 'node2'], ['4M', '7M'], tags={}, labels={'label1': 'test'})

    result = mod.serialize_gfa(gfa_graph)
    expected = '{"paths": {"000000F": {"labels": {"label1": "test"}, "nodes": ["node1", "node2"], "tags": {}, "name": "000000F", "cigars": ["4M", "7M"]}}, "nodes": {"node1": {"labels": {}, "seq": "ACTGAAA", "name": "node1", "len": 7, "tags": {}}, "node2": {"labels": {}, "seq": "AAACCCGGGT", "name": "node2", "len": 10, "tags": {}}}, "edges": {"(\'node1\', \'node2\')": {"labels": {}, "v_orient": "+", "tags": {}, "v_start": 4, "cigar": "*", "w_end": 3, "w_start": 0, "w_orient": "+", "name": "edge1", "v_end": 7, "w": "node2", "v": "node1"}}}'

    assert(result == expected)

def test_deserialize():
    gfa_graph = mod.GFAGraph()
    gfa_graph.add_node('node1', 7, 'ACTGAAA', tags={}, labels={})
    gfa_graph.add_node('node2', 10, 'AAACCCGGGT', tags={}, labels={})
    gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 7, 0, 3, '*', tags={}, labels={})
    gfa_graph.add_path('000000F', ['node1', 'node2'], ['4M', '7M'], tags={}, labels={'label1': 'test'})

    dump = mod.serialize_gfa(gfa_graph)

    fp_in = StringIO(dump)

    result = mod.deserialize_gfa(fp_in)

    assert (result.nodes == gfa_graph.nodes)
    assert (result.edges == gfa_graph.edges)
    assert (result.paths == gfa_graph.paths)
