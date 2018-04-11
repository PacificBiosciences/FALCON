import argparse
import os
import sys
import json

from falcon_kit.gfa_graph import GFAGraph, serialize_gfa, deserialize_gfa
import falcon_kit.mains.collect_pread_gfa
import falcon_kit.tiling_path

def run(fp_out, p_ctg_tiling_path, a_ctg_tiling_path,
        p_ctg_fasta, a_ctg_fasta,
        write_contigs,
        min_p_len, min_a_len, only_these_contigs):

    gfa_graph = GFAGraph()

    # Load the primary and associate contig files.
    p_ctg_dict = falcon_kit.mains.collect_pread_gfa.load_seqs(p_ctg_fasta, (not write_contigs))
    p_ctg_lens = {key: val[0] for key, val in p_ctg_dict.iteritems()}
    p_ctg_seqs = {key: val[1] for key, val in p_ctg_dict.iteritems()}

    a_ctg_dict = falcon_kit.mains.collect_pread_gfa.load_seqs(a_ctg_fasta, (not write_contigs))
    a_ctg_lens = {key: val[0] for key, val in a_ctg_dict.iteritems()}
    a_ctg_seqs = {key: val[1] for key, val in a_ctg_dict.iteritems()}

    # Create whitelists for filtering contigs.
    p_ctg_whitelist = set(p_ctg_seqs.keys())
    a_ctg_whitelist = set([key for key in a_ctg_seqs.keys()])
    if only_these_contigs:
        p_ctg_whitelist = set(open(only_these_contigs).read().splitlines()) & set(p_ctg_whitelist)
        a_ctg_whitelist = set([key for key in a_ctg_seqs.keys() if key.split('-')[0].split('_')[0] in p_ctg_whitelist])

    # Load the tiling paths and assign coordinates.
    p_paths = falcon_kit.tiling_path.load_tiling_paths(p_ctg_tiling_path, whitelist_seqs=p_ctg_whitelist, contig_lens=p_ctg_lens)
    a_paths = falcon_kit.tiling_path.load_tiling_paths(a_ctg_tiling_path, whitelist_seqs=a_ctg_whitelist, contig_lens=a_ctg_lens)

    # Find the associate contig placement. `a_placement` is a dict:
    #   placement[p_ctg_id][a_ctg_id] = (start, end, p_ctg_id, a_ctg_id, first_node, last_node)
    a_placement = falcon_kit.tiling_path.find_a_ctg_placement(p_paths, a_paths)

    # Add the nodes.
    for ctg_id, tiling_path in p_paths.iteritems():
        gfa_graph.add_node(ctg_id, p_ctg_lens[ctg_id], p_ctg_seqs[ctg_id])
    for ctg_id, tiling_path in a_paths.iteritems():
        gfa_graph.add_node(ctg_id, a_ctg_lens[ctg_id], a_ctg_seqs[ctg_id])

    for p_ctg_id, a_dict in a_placement.iteritems():
        for a_ctg_id, placement in a_dict.iteritems():
            start, end, p_ctg_id, a_ctg_id, first_node, last_node = placement

            a_ctg_len = a_ctg_lens[a_ctg_id]

            # edge_name = 'edge-%d-out-%s-to-%s' % (len(gfa_graph.edges), a_ctg_id, p_ctg_id)
            edge_name = 'edge-%d' % (len(gfa_graph.edges))
            gfa_graph.add_edge(edge_name, p_ctg_id, '+', a_ctg_id, '+', start, start, 0, 0, '*', tags = {}, labels = {})

            # edge_name = 'edge-%d-in-%s-to-%s' % (len(gfa_graph.edges), a_ctg_id, p_ctg_id)
            edge_name = 'edge-%d' % (len(gfa_graph.edges))
            gfa_graph.add_edge(edge_name, a_ctg_id, '+', p_ctg_id, '+', a_ctg_len, a_ctg_len, end, end, '*', tags = {}, labels = {})

    fp_out.write(serialize_gfa(gfa_graph))
    fp_out.write('\n')

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Generates GFA output (on stdout) from FALCON's assembly.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--p-ctg-tiling-path', type=str, default='p_ctg_tiling_path',
                        help='location of the p_ctg tiling path file')
    parser.add_argument('--a-ctg-tiling-path', type=str, default='a_ctg_tiling_path',
                        help='location of the a_ctg tiling path file')
    parser.add_argument('--p-ctg-fasta', type=str, default='p_ctg.fa',
                        help='path to the primary contigs file')
    parser.add_argument('--a-ctg-fasta', type=str, default='a_ctg.fa',
                        help='path to the associate contigs file')
    parser.add_argument('--write-contigs', '-c', action='store_true',
                        help="output contig sequences as S lines")
    parser.add_argument('--min-p-len', type=int, default=0,
                        help='primary contig paths with length smaller than this will not be reported')
    parser.add_argument('--min-a-len', type=int, default=0,
                        help='associate contig paths with length smaller than this will not be reported')
    parser.add_argument('--only-these-contigs', type=str, default='',
                        help='limit output to specified contigs listed in file (one per line)')
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)

    run(sys.stdout, **vars(args))

if __name__ == '__main__':  # pragma: no cover
    main()
