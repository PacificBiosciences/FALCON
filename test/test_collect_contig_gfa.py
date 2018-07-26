import falcon_kit.mains.collect_contig_gfa as mod
import helpers
import pytest
import os
import falcon_kit.tiling_path
from falcon_kit.FastaReader import FastaReader
from test_tiling_path import test_placement_1_p_path_as_text, test_placement_1_a_path_as_text
import json

import random
random.seed(1234567)

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

# def load_seq(fasta_fn):
#     return {r.name.split()[0]: r.sequence.upper() for r in FastaReader(fasta_fn)}

# def path_to_seq(preads, path_as_split_lines, with_first_read):
#     ret = ''

#     if len(path_as_split_lines) == 0:
#         return ret

#     if with_first_read:
#         ctg_id, v, w, wrid, sp, tp = path_as_split_lines[0][0:6]
#         vrid, vorient = v.split(':')
#         ret += preads[vrid] if vorient == 'E' else "".join([RCMAP[c] for c in preads[vrid][::-1]])

#     for edge in path_as_split_lines:
#         ctg_id, v, w, wrid, sp, tp = edge[0:6]
#         sp, tp = int(sp), int(tp)
#         ret += preads[wrid][sp:tp] if sp < tp else "".join([RCMAP[c] for c in preads[wrid][sp:tp:-1]])

#     ret = ''.join([str(val) for val in ret])

#     return ret

def generate_seq(seq_len):
    return ''.join([random.choice(['A', 'C', 'T', 'G']) for i in xrange(seq_len)])

def setup_test(fp_p_ctg, fp_p_tp, fp_a_ctg, fp_a_tp):
    """
    Configures test files with predefined values.
    This includes p and a tiling paths and contig sequences, and also returns
    the values as dicts.
    Expected value is returned, for the case where sequences should be
    omitted.
    """

    p_tp_as_text = test_placement_1_p_path_as_text
    a_tp_as_text = test_placement_1_a_path_as_text

    all_seqs = {}

    # Write sample p_ctg_tiling_path.
    fp_p_tp.write(p_tp_as_text)

    # Load and parse the tiling path. Validity is tested in test_tiling_path.
    p_ctg_tps = falcon_kit.tiling_path.load_tiling_paths(str(fp_p_tp), contig_lens=None, whitelist_seqs=None)

    # Generate random p_ctg sequences write to file.
    p_ctg_fasta_text = ''
    for key, tp in p_ctg_tps.iteritems():
        all_seqs[key] = generate_seq(tp.contig_len)
        p_ctg_fasta_text += '>%s\n%s\n' % (key, all_seqs[key])
    fp_p_ctg.write(p_ctg_fasta_text)


    # Write sample p_ctg_tiling_path.
    fp_a_tp.write(a_tp_as_text)

    # Load and parse the tiling path. Validity is tested in test_tiling_path.
    a_ctg_tps = falcon_kit.tiling_path.load_tiling_paths(str(fp_a_tp), contig_lens=None, whitelist_seqs=None)

    # Generate random p_ctg sequences write to file.
    a_ctg_fasta_text = ''
    for key, tp in a_ctg_tps.iteritems():
        all_seqs[key] = generate_seq(tp.contig_len)
        a_ctg_fasta_text += '>%s\n%s\n' % (key, all_seqs[key])
    fp_a_ctg.write(a_ctg_fasta_text)

    expected = \
            {
                "nodes": {
                    "000000F": {
                        "labels": {}, "seq": "*", "name": "000000F", "len": 71411, "tags": {}
                        },
                    "000001F": {
                        "labels": {}, "seq": "*", "name": "000001F", "len": 33726, "tags": {}
                        },
                    "000002F": {
                        "labels": {}, "seq": "*", "name": "000002F", "len": 8473, "tags": {}
                        },
                    "000000F-1": {
                        "labels": {}, "seq": "*", "name": "000000F-1", "len": 29405, "tags": {}
                        }
                },
                "edges": {
                    "('000000F', '000000F-1')": {
                        "labels": {},
                        "tags": {},
                        "v": "000000F",
                        "v_orient": "+",
                        "w": "000000F-1",
                        "w_orient": "+",
                        "v_start": 43849,
                        "v_end": 43849,
                        "w_start": 0,
                        "w_end": 0,
                        "cigar": "*",
                        "name": "edge-0"},
                    "('000000F-1', '000000F')": {
                        "labels": {},
                        "tags": {},
                        "v": "000000F-1",
                        "v_orient": "+",
                        "w": "000000F",
                        "w_orient": "+",
                        "v_start": 29405,
                        "v_end": 29405,
                        "w_start": 67383,
                        "w_end": 67383,
                        "cigar": "*",
                        "name": "edge-1"},
                },
                "paths": {},
            }

    return all_seqs, p_ctg_tps, a_ctg_tps, expected

def test_main_1(tmpdir, capsys):
    """
    Test without writing contig sequences.
    """

    p_ctg_file = tmpdir.join('p_ctg.fa')
    p_tp_file = tmpdir.join('p_ctg_tiling_path')
    a_tp_file = tmpdir.join('a_ctg_tiling_path')
    a_ctg_file = tmpdir.join('a_ctg.fa')
    all_seqs, p_ctg_tps, a_ctg_tps, expected = setup_test(p_ctg_file, p_tp_file, a_ctg_file, a_tp_file)

    argv = ['prog',
            '--p-ctg-tiling-path', str(p_tp_file),
            '--a-ctg-tiling-path', str(a_tp_file),
            '--p-ctg-fasta', str(p_ctg_file),
            '--a-ctg-fasta', str(a_ctg_file),
            # '--write-contigs',
            '--min-p-len', '0',
            '--min-a-len', '0',
            # '--only-these-contigs'
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    result = json.loads(out)

    assert(result == expected)

def test_main_2(tmpdir, capsys):
    """
    Test output with sequences.
    """

    p_ctg_file = tmpdir.join('p_ctg.fa')
    p_tp_file = tmpdir.join('p_ctg_tiling_path')
    a_tp_file = tmpdir.join('a_ctg_tiling_path')
    a_ctg_file = tmpdir.join('a_ctg.fa')
    all_seqs, p_ctg_tps, a_ctg_tps, expected = setup_test(p_ctg_file, p_tp_file, a_ctg_file, a_tp_file)

    # This test expects that the nodes are initialized with sequences, so
    # set the reference to the generated sequences here.
    for key, node in expected['nodes'].iteritems():
        node['seq'] = all_seqs[key]

    argv = ['prog',
            '--p-ctg-tiling-path', str(p_tp_file),
            '--a-ctg-tiling-path', str(a_tp_file),
            '--p-ctg-fasta', str(p_ctg_file),
            '--a-ctg-fasta', str(a_ctg_file),
            '--write-contigs',
            '--min-p-len', '0',
            '--min-a-len', '0',
            # '--only-these-contigs'
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    result = json.loads(out)

    assert(result == expected)

def test_main_3(tmpdir, capsys):
    """
    Test the whitelist.
    """

    p_ctg_file = tmpdir.join('p_ctg.fa')
    p_tp_file = tmpdir.join('p_ctg_tiling_path')
    a_tp_file = tmpdir.join('a_ctg_tiling_path')
    a_ctg_file = tmpdir.join('a_ctg.fa')
    all_seqs, p_ctg_tps, a_ctg_tps, expected = setup_test(p_ctg_file, p_tp_file, a_ctg_file, a_tp_file)

    only_these_contigs = set(['000002F'])
    whitelist_file = tmpdir.join('only_these_contigs')
    whitelist_file.write('\n'.join(only_these_contigs) + '\n')

    # Adjust the expected results.
    # Remove any node not in the whitelist.
    blacklist_nodes = [key for key in (p_ctg_tps.keys() + a_ctg_tps.keys()) if key not in only_these_contigs]
    for key in blacklist_nodes:
        expected['nodes'].pop(key, None)
    # Remove any edge not in the whitelist.
    blacklist_edges = [key for key, edge in expected['edges'].iteritems() if edge['v'] not in only_these_contigs or edge['w'] not in only_these_contigs]
    for key in blacklist_edges:
            expected['edges'].pop(key, None)

    argv = ['prog',
            '--p-ctg-tiling-path', str(p_tp_file),
            '--a-ctg-tiling-path', str(a_tp_file),
            '--p-ctg-fasta', str(p_ctg_file),
            '--a-ctg-fasta', str(a_ctg_file),
            # '--write-contigs',
            '--min-p-len', '0',
            '--min-a-len', '0',
            '--only-these-contigs', str(whitelist_file),
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    result = json.loads(out)

    assert(result == expected)

def test_main_4(tmpdir, capsys):
    """
    Test for an empty input.
    """

    p_ctg_file = tmpdir.join('p_ctg.fa')
    p_tp_file = tmpdir.join('p_ctg_tiling_path')
    a_tp_file = tmpdir.join('a_ctg_tiling_path')
    a_ctg_file = tmpdir.join('a_ctg.fa')
    p_ctg_file.write('')
    p_tp_file.write('')
    a_tp_file.write('')
    a_ctg_file.write('')

    expected = {'nodes': {}, 'edges': {}, 'paths': {}}

    argv = ['prog',
            '--p-ctg-tiling-path', str(p_tp_file),
            '--a-ctg-tiling-path', str(a_tp_file),
            '--p-ctg-fasta', str(p_ctg_file),
            '--a-ctg-fasta', str(a_ctg_file),
            # '--write-contigs',
            '--min-p-len', '0',
            '--min-a-len', '0',
            # '--only-these-contigs'
            ]
    mod.main(argv)
    out, err = capsys.readouterr()

    result = json.loads(out)

    assert(result == expected)
