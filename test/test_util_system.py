import falcon_kit.util.system as mod
import helpers
import pytest
import os

def test_find_files():
    # Testing a normal case.
    root_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')
    result = list(mod.find_files(root_dir, '*.gexf'))
    expected = [os.path.join(root_dir, 'sg.gexf')]
    assert(result == expected)

    # Testing a normal case from one level up.
    root_dir = helpers.get_test_data_dir()
    result = list(mod.find_files(root_dir, '*.gexf'))
    expected = [os.path.join(root_dir, 'gfa-1', 'sg.gexf')]
    assert(result == expected)

    # Testing an empty pattern.
    root_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')
    result = list(mod.find_files(root_dir, ''))
    expected = []
    assert(result == expected)

    # Testing the lookup of everything.
    root_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')
    result = mod.find_files(root_dir, '*')
    expected = ['a_ctg.fa', 'a_ctg_tiling_path', 'p_ctg.fa', 'p_ctg_tiling_path', 'ctg_paths',
                'preads4falcon.fasta', 'sg.gexf', 'sg_edges_list', 'utg_data',
                'expected-1-sg-r-c.gfa', 'expected-2-tiling-r-c.gfa', 'expected-3-tiling-no_r-c.gfa',
                'expected-4-tiling-no_r-no_c.gfa', 'expected-5-sg-no_r-no_c.gfa',
                'expected-6-tiling-no_r-no_c-minlen.gfa', 'expected-7-nx-no_r-no_c.gfa',
                'expected-8-nx-tiling-no_r-no_c.gfa', 'expected-9-nx-tiling-r-c.gfa',
                ]
    expected = [os.path.join(root_dir, val) for val in expected]
    assert(sorted(result) == sorted(expected))

    # Testing an empty directory and an empty pattern.
    result = list(mod.find_files('', ''))
    expected = []
    assert(result == expected)

    # Testing a non existent pattern.
    root_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')
    result = list(mod.find_files(root_dir, 'Wubalubadubdub'))
    expected = []
    assert(result == expected)
