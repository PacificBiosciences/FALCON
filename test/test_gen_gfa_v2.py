import falcon_kit.mains.gen_gfa_v2 as mod
import helpers
import pytest
import os

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

def test_main_1(tmpdir, capsys):
    # test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    gfa_graph = mod.GFAGraph()
    gfa_graph.add_node('node1', 7, 'ACTGAAA', tags={}, labels={})
    gfa_graph.add_node('node2', 10, 'AAACCCGGGT', tags={}, labels={})
    gfa_graph.add_edge('edge1', 'node1', '+', 'node2', '+', 4, 7, 0, 3, '*', tags={}, labels={})
    gfa_graph.add_path('000000F', ['node1', 'node2'], ['4M', '7M'], tags={}, labels={})

    gfa_json_file = tmpdir.join('graph.gfa.json')
    gfa_json_file.write(mod.serialize_gfa(gfa_graph))

    argv = ['prog',
            str(gfa_json_file),
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    expected = """H	VN:Z:2.0
S	node1	7	ACTGAAA
S	node2	10	AAACCCGGGT
E	edge1	node1+	node2+	4	7$	0	3	*
"""

    # Custom tags can be in an arbitrary order and hard to compare.
    # They might also change between versions, so just remove
    # them manually here and compare only the important fields.
    result = []
    for line in out.splitlines():
        line = line.strip()
        sl = line.split()
        if len(line) == 0:
            continue
        if line[0] == 'H':
            result.append('\t'.join(sl[0:2]))
        elif line[0] == 'S':
            result.append('\t'.join(sl[0:4]))
        elif line[0] == 'E':
            result.append('\t'.join(sl[0:9]))

    result_str = '\n'.join(result) + '\n'

    assert(result_str == expected)
