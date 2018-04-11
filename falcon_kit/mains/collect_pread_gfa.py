import argparse
import os
import sys
import json

from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader
from falcon_kit.gfa_graph import *
import falcon_kit.tiling_path

def load_seqs(fasta_fn, store_only_seq_len):
    """
    If store_only_seq_len is True, then the seq is discarded and
    only it's length stored.
    """
    seqs = {}
    f = FastaReader(fasta_fn)
    if store_only_seq_len == False:
        for r in f:
            seqs[r.name.split()[0]] = (len(r.sequence), r.sequence.upper())
    else:
        for r in f:
            seqs[r.name.split()[0]] = (len(r.sequence), '*')
    return seqs

def load_pread_overlaps(fp_in):
    preads_overlap_dict = {}
    for line in fp_in:
        sl = line.strip().split()
        if len(sl) < 13:
            continue
        # Example line: 000000009 000000082 -3004 99.90 0 4038 7043 7043 1 6488 9492 9492 overlap 000000F.5000003.0 000000F.5000003.0
        preads_overlap_dict[(sl[0], sl[1])] = sl[0:4] + [int(val) for val in sl[4:12]] + sl[12:]
        
        # Overlaps are not always symmetrically represented in the preads.ovl for some reason, so add the
        # reverse overlap here as well, but do not overwrite existing (just to be safe).
        if (sl[1], sl[0]) not in preads_overlap_dict:
            preads_overlap_dict[(sl[1], sl[0])] = [sl[1], sl[0], sl[2], sl[3]] + [int(val) for val in sl[8:12]] + [int(val) for val in sl[4:8]] + sl[12:]

    return preads_overlap_dict

def load_sg_edges(fp_in):
    """
    Loads all sg_edges_list so that haplotig paths can be reversed if needed.
    with open(os.path.join(fc_asm_path, "sg_edges_list"), 'r') as fp:
        sg_edges_dict = load_sg_edges(fp)
    """
    sg_edges_dict = {}
    for line in fp_in:
        sl = line.strip().split()
        if len(sl) < 8:
            continue
        # Example line: 000000512:B 000000679:E 000000679 4290 7984 4290 99.95 TR
        sg_edges_dict[(sl[0], sl[1])] = sl[0:3] + [int(val) for val in sl[3:6]] + [float(sl[6])] + sl[7:]
    return sg_edges_dict

def add_node(gfa_graph, v, preads_dict):
    v_name, v_orient = v.split(':')
    v_len, v_seq = preads_dict[v_name]
    gfa_graph.add_node(v_name, v_len, v_seq)

def add_edge(gfa_graph, v, w, edge_split_line, preads_overlap_dict, sg_edges_dict):
    edge_name = 'edge-%d' % (len(gfa_graph.edges))
    v_name, v_orient = v.split(':')
    w_name, w_orient = w.split(':')
    v_orient = '+' if v_orient == 'E' else '-'
    w_orient = '+' if w_orient == 'E' else '-'
    cigar = '*'

    # Get the SG edge and the overlap, and set the tags and labels.
    sg_edge = sg_edges_dict[(v, w)]
    overlap = preads_overlap_dict[(v_name, w_name)]
    labels = {'tp': edge_split_line, 'sg_edge': sg_edge, 'overlap': overlap}
    tags = {}

    # Example overlap:
    #   000000001 000000170 -6104 99.75 0 1909 8010 8010 1 1250 7354 7354 overlap 000000F.5000003.0 000000F.5000003.0
    # Handle the overlap coordinates - GFA format requires the coordinates to be with
    # respect to the fwd strand, and the M4 format reports overlaps on the
    # strand of the alignment.
    _, _, score, idt, v_rev, v_start, v_end, v_len, w_rev, w_start, w_end, w_len = overlap[0:12]
    if v_rev == 1:
        v_start, v_end = v_end, v_start
        v_start = v_len - v_start
        v_end = v_len - v_end
    if w_rev == 1:
        w_start, w_end = w_end, w_start
        w_start = w_len - w_start
        w_end = w_len - w_end

    gfa_graph.add_edge(edge_name, v_name, v_orient, w_name, w_orient, v_start, v_end, w_start, w_end, cigar, tags = tags, labels = labels)

def add_tiling_paths_to_gfa(gfa_graph, tiling_paths, preads_dict, preads_overlap_dict, sg_edges_dict):
    # Add nodes.
    for ctg_id, tiling_path in tiling_paths.iteritems():
        for edge in tiling_path.edges:
            add_node(gfa_graph, edge.v, preads_dict)
            add_node(gfa_graph, edge.w, preads_dict)

    # Add edges.
    for ctg_id, tiling_path in tiling_paths.iteritems():
        for edge in tiling_path.edges:
            add_edge(gfa_graph, edge.v, edge.w, edge.get_split_line(), preads_overlap_dict, sg_edges_dict)

    # Add path.
    for ctg_id, tiling_path in tiling_paths.iteritems():
        path_nodes = []
        path_cigars = []
        if len(tiling_path.edges) == 0:
            continue

        # Add the first node to the path.
        v = tiling_path.edges[0].v
        v_name, v_orient = v.split(':')
        cigar = '%dM' % (tiling_path.coords[v]) # This will be 0 if the contig is improper, and length of v otherwise.
        path_nodes.append(v_name)
        path_cigars.append(cigar)

        # Add the rest of the nodes.
        for edge in tiling_path.edges:
            w_name, w_orient = edge.w.split(':')
            cigar = '%dM' % (abs(edge.e - edge.b))
            path_nodes.append(w_name)
            path_cigars.append(cigar)

        gfa_graph.add_path(ctg_id, path_nodes, path_cigars)

def add_string_graph_to_gfa(gfa_graph, sg_edges_list, utg_data, ctg_paths, preads_dict, preads_overlap_dict, sg_edges_dict):
    asm_graph = AsmGraph(sg_edges_list, utg_data, ctg_paths)

    for v, w in asm_graph.sg_edges:
        add_node(gfa_graph, v, preads_dict)
        add_node(gfa_graph, w, preads_dict)

    for v, w in asm_graph.sg_edges:
        edge_data = asm_graph.sg_edges[(v, w)]
        if edge_data[-1] != 'G':
            continue
        add_edge(gfa_graph, v, w, edge_data, preads_overlap_dict, sg_edges_dict)

def run(fp_out, p_ctg_tiling_path, a_ctg_tiling_path,
                      preads_fasta, p_ctg_fasta, a_ctg_fasta,
                      sg_edges_list, preads_ovl, utg_data, ctg_paths,
                      add_string_graph, write_reads,
                      min_p_len, min_a_len, only_these_contigs):
    """
    This method produces a GFAGraph object containing info required
    to write both the GFA-1 and GFA-2 formatted assemblies.
    However, it does not write the GFA formats directly, but instead
    dumps a JSON file to disk.
    The JSON file is converted to a GFA-1 or a GFA-2 with outside scripts.

    The graphical output is produced from either the entire string
    graph (only the non-filtered edges are considered) or from only
    the tiling paths. String graph can show the neighborhood of contig
    breaks, whereas the tiling path output is more sparse.
    Output is written to stdout.
    """

    gfa_graph = GFAGraph()

    # Load preads.
    preads_dict = load_seqs(preads_fasta, (not write_reads))

    # Load the pread overlaps
    with open(preads_ovl, 'r') as fp:
        preads_overlap_dict = load_pread_overlaps(fp)

    # Load the SG edges.
    with open(sg_edges_list, 'r') as fp:
        sg_edges_dict = load_sg_edges(fp)

    # Load the primary and associate contig files.
    p_ctg_seqs = load_seqs(p_ctg_fasta, True)
    a_ctg_seqs = load_seqs(a_ctg_fasta, True)

    # Collect the sequence lengths from the above dicts.
    p_ctg_lens = {key: val[0] for key, val in p_ctg_seqs.iteritems()}
    a_ctg_lens = {key: val[0] for key, val in a_ctg_seqs.iteritems()}

    # Create whitelists for filtering contigs.
    p_ctg_whitelist = set(p_ctg_seqs.keys())
    a_ctg_whitelist = set([key for key in a_ctg_seqs.keys()])
    if only_these_contigs:
        p_ctg_whitelist = set(open(only_these_contigs).read().splitlines()) & set(p_ctg_whitelist)
        a_ctg_whitelist = set([key for key in a_ctg_seqs.keys() if key.split('-')[0].split('_')[0] in p_ctg_whitelist])

    # Load the tiling paths and assign coordinates.
    p_paths = falcon_kit.tiling_path.load_tiling_paths(p_ctg_tiling_path, whitelist_seqs=p_ctg_whitelist, contig_lens=p_ctg_lens)
    a_paths = falcon_kit.tiling_path.load_tiling_paths(a_ctg_tiling_path, whitelist_seqs=a_ctg_whitelist, contig_lens=a_ctg_lens)

    add_tiling_paths_to_gfa(gfa_graph, p_paths, preads_dict, preads_overlap_dict, sg_edges_dict)
    add_tiling_paths_to_gfa(gfa_graph, a_paths, preads_dict, preads_overlap_dict, sg_edges_dict)

    if add_string_graph:
        add_string_graph_to_gfa(gfa_graph, sg_edges_list, utg_data, ctg_paths, preads_dict, preads_overlap_dict, sg_edges_dict)

    fp_out.write(serialize_gfa(gfa_graph))
    fp_out.write('\n')

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Generates GFA output (on stdout) from FALCON's assembly.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--p-ctg-tiling-path', type=str, default='p_ctg_tiling_path',
                        help='location of the p_ctg tiling path file')
    parser.add_argument('--a-ctg-tiling-path', type=str, default='a_ctg_tiling_path',
                        help='location of the a_ctg tiling path file')
    parser.add_argument('--preads-fasta', type=str, default='preads4falcon.fasta',
                        help='path to the preads4falcon.fasta file')
    parser.add_argument('--p-ctg-fasta', type=str, default='p_ctg.fa',
                        help='path to the primary contigs file')
    parser.add_argument('--a-ctg-fasta', type=str, default='a_ctg.fa',
                        help='path to the associate contigs file')
    parser.add_argument('--sg-edges-list', type=str, default='sg_edges_list',
                        help='string graph edges file from Falcon assembly')
    parser.add_argument('--preads-ovl', type=str, default='preads.ovl',
                        help='the preads overlap file')
    parser.add_argument('--utg-data', type=str,
                        default='utg_data', help='unitig data file from Falcon')
    parser.add_argument('--ctg-paths', type=str, default='ctg_paths',
                        help='contig paths file from Falcon assembly')
    parser.add_argument('--add-string-graph', action='store_true',
                        help="in addition to tiling paths, output other edges and nodes from the final string graph")
    parser.add_argument('--write-reads', '-r', action='store_true',
                        help="output read sequences in S lines")
    # parser.add_argument('--write-contigs', '-c', action='store_true',
    #                     help="output contig sequences as S lines")
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
