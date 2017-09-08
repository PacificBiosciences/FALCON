import argparse
import logging
import os
import sys

GFA_H_TAG = 'H'
GFA_S_TAG = 'S'
GFA_L_TAG = 'L'
GFA_P_TAG = 'P'
GFA_ORIENT_FWD = '+'
GFA_ORIENT_REV = '-'
GFA_SEQ_UNKNOWN = '*'
GFA_LINK_CIGAR_UNKNOWN = '*'

from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader

def calc_node_coords(tiling_path):
    """
    For a single tiling path (tiling_path is a list
    of edges for a particular contig) calculates the
    genomic coordinate of every node in the path.
    """
    if not tiling_path:
        return {}, 0
    coord_map = {}
    contig_len = 0
    v0, w0, b, e, l, idt, etype = tiling_path[0]
    coord_map[v0] = 0
    for v, w, b, e, l, idt, etype in tiling_path:
        if v not in coord_map:
            raise Exception('Tiling path is not in sorted order. Node "{v}" does not yet have an assigned coordinate.'.format(v=v))
        if w in coord_map and w != v0:
            raise Exception('Tiling path is not linear. Node "{w}" has multiple branches.'.format(w=w))
        coord = coord_map[v]
        coord += abs(int(b) - int(e))
        coord_map[w] = coord
        contig_len = max(contig_len, coord)
    return coord_map, contig_len

def load_tiling_paths_from_stream(fp, type_):
    """
    Loads the tiling paths from a specified tiling path file stream.
    Parameter type_ is intended to mark the type of the tiling path
    file (e.g. 'P' for primary contigs and 'A' for associate contigs).
    """
    paths = {}
    edge_to_ctg = {}
    if not fp:
        return paths, edge_to_ctg
    for row in fp:
        row = row.strip().split()
        ctg_id, v, w, edge_rid, b, e, l, idt  = row[:9]
        paths.setdefault(ctg_id, [])
        paths[ctg_id].append( (v, w, int(b), int(e), int(l), float(idt), type_) )
        ctg_id = ctg_id.split('-')[0] #get the primary contig id
        edge_to_ctg[(v, w)] = ctg_id, type_
    return paths, edge_to_ctg

def load_tiling_paths(tiling_path_file, type_):
    """
    Loads the tiling paths from a specified tiling path file.
    Parameter type_ is intended to mark the type of the tiling path
    file (e.g. 'P' for primary contigs and 'A' for associate contigs).
    """
    paths = {}
    edge_to_ctg = {}
    with open(tiling_path_file) as f:
        paths, edge_to_ctg = load_tiling_paths_from_stream(f, type_)
    return paths, edge_to_ctg

def calc_tiling_paths_len(tiling_paths):
    """
    Calculates the length of each tiling path in the given
    tiling_paths dict. Inherently, it computes the starting position
    of each node in the tiling path as well.
    """
    path_coords = {}
    paths_len = {}
    for ctg_id, edges in tiling_paths.iteritems():
        path_coords[ctg_id], paths_len[ctg_id] = calc_node_coords(edges)
    return path_coords, paths_len

def filter_tiling_paths_by_len(tiling_paths, paths_len, min_len):
    """
    Given a dict of tiling_paths (key is contig name, value is a
    list of tiling path edges) creates and returns a new dict with
    tiling paths which are larger than the specified minimum length.
    """
    ret_paths = {}
    for ctg_id, edges in tiling_paths.iteritems():
        if ctg_id not in paths_len: continue
        plen = paths_len[ctg_id]
        if plen >= min_len:
            ret_paths[ctg_id] = edges
    return ret_paths

def get_gfa_links_from_sg(G_asm, edge_to_ctg):
    """
    Iterates through string graph edges and selects only
    those edges which made it to the final assembly graph
    (these edges are marked as 'G').
    This method creates a list of links for the GFA output.
    """
    read_in_graph = set()
    read_pairs = set()
    link_lines = []
    for v, w in G_asm.sg_edges:
        edge_data = G_asm.sg_edges[(v, w)]
        """000084959:E 000376804:B 000376804 25867 0 1766 99.55 TR"""
        """edge_data looks like this: (('000024576', 27607, 0), 6377, 99.66, 'G')"""
        if edge_data[-1] == 'G':
            r1, r1end = v.split(':')
            r2, r2end = w.split(':')
            rp = [r1, r2]
            rp.sort()
            rp = tuple(rp)
            if rp in read_pairs:  #GFA model read rather than the end, so we need to collapse dual edge
                continue

            read_pairs.add(rp)
            read_in_graph.add(r1)
            read_in_graph.add(r2)
            if r1end == 'E':
                o1 = GFA_ORIENT_FWD
            else:
                o1 = GFA_ORIENT_REV
            if r2end == 'E':
                o2 = GFA_ORIENT_FWD
            else:
                o2 = GFA_ORIENT_REV
            overlap_begin = int(edge_data[0][1])
            overlap_end = int(edge_data[0][2])
            overlap_length = int(edge_data[1])
            overlap_idt = float(edge_data[2])
            ctg_id = edge_to_ctg.get((v, w), ('NA', 'NA'))
            link_lines.append(  '\t'.join([GFA_L_TAG, r1, o1, r2, o2, GFA_LINK_CIGAR_UNKNOWN,
                                'ol:i:%d' % overlap_length, 'oi:f:%.1f' % overlap_idt,
                                'ob:i:%d' % overlap_begin, 'oe:i:%d' % overlap_end,
                                'ci:A:%s-%s' % ctg_id]) )

    return read_in_graph, link_lines

def get_gfa_links_from_tiling_paths(tiling_paths, edge_to_ctg):
    """
    This method creates a list of links for the GFA output based on
    the tiling paths loaded from p_ctg_tiling_path and a_ctg_tiling_path.
    """
    read_in_graph = set()
    read_pairs = set()
    link_lines = []
    for ctg_header, tiling_path in tiling_paths.iteritems():
        for v, w, b, e, l, idt, etype in tiling_path:
            r1, r1end = v.split(':')
            r2, r2end = w.split(':')
            rp = [r1, r2]
            rp.sort()
            rp = tuple(rp)
            if rp in read_pairs:  #GFA model read rather than the end, so we need to collapse dual edge
                continue

            read_pairs.add(rp)
            read_in_graph.add(r1)
            read_in_graph.add(r2)
            if r1end == 'E':
                o1 = GFA_ORIENT_FWD
            else:
                o1 = GFA_ORIENT_REV
            if r2end == 'E':
                o2 = GFA_ORIENT_FWD
            else:
                o2 = GFA_ORIENT_REV
            ctg_id = edge_to_ctg.get((v, w), ('NA', 'NA'))
            link_lines.append(  '\t'.join([GFA_L_TAG, r1, o1, r2, o2, GFA_LINK_CIGAR_UNKNOWN,
                                'ol:i:%d' % l, 'oi:f:%.1f' % idt,
                                'ob:i:%d' % b, 'oe:i:%d' % e,
                                'ci:A:%s-%s' % ctg_id]) )

    return read_in_graph, link_lines

def format_gfa_path_line(ctg_name, tiling_path, seq_len):
    """
    Takes a tiling path and formulates a P line for GFA.
    """
    if not tiling_path:
        return []
    v, w, b, e, l, idt, etype = tiling_path[0]
    rn, rend = v.split(':')
    o = GFA_ORIENT_FWD if rend == 'E' else GFA_ORIENT_REV
    segs = [rn+o]
    segs_cigar = ['%dM' % (seq_len[rn])]        # This is not part of p_ctg.fa.
    for v, w, b, e, l, idt, etype in tiling_path:
        rn, rend = w.split(':')
        o = GFA_ORIENT_FWD if rend == 'E' else GFA_ORIENT_REV
        segs.append(rn + o)
        l = abs(b - e)
        segs_cigar.append('%dM' % (l))
    out = [GFA_P_TAG, ctg_name, ','.join(segs), ','.join(segs_cigar)]
    return out;

def load_seqs(fasta_file, seqs_to_load, load_bases):
    """
    Opens a FASTA file and loads sequences with headers
    present in the filter set (seqs_to_load).
    Returns the sequences (if load_bases == True) and their
    lengths (returned always).
    If seqs_to_load is an empty set, all sequences
    will be loaded.
    """
    f = FastaReader(fasta_file)
    seq_len = {}
    seq_bases = {}
    for r in f:
        if seqs_to_load and r.name not in seqs_to_load:
            continue
        seq_len[r.name] = len(r.sequence)
        if (load_bases == True):
            seq_bases[r.name] = r.sequence
    return seq_bases, seq_len

def gfa_from_assembly(out, p_ctg_tiling_path, a_ctg_tiling_path, \
                      preads_fasta, p_ctg_fasta, a_ctg_fasta, \
                      sg_edges_list, utg_data, ctg_paths, \
                      tiling, write_reads, write_contigs, \
                      min_p_len, min_a_len):
    """
    This method produces the GFA-1 formatted output of the
    FALCON assembly.
    The graphical output is produced from either the entire string
    graph (only the non-filtered edges are considered) or from only
    the tiling paths. String graph can show the neighborhood of contig
    breaks, whereas the tiling path output is more sparse.
    Output is written to stdout.
    """

    # Initialize the output stream.
    fp_out = sys.stdout if out == '-' else open(out, 'w')

    # Load and filter primary contig paths.
    p_path, p_edge_to_ctg = load_tiling_paths(p_ctg_tiling_path, 'P');
    _, p_ctg_len = calc_tiling_paths_len(p_path)
    p_path = filter_tiling_paths_by_len(p_path, p_ctg_len, min_p_len)

    # Load and filter associate contig paths.
    a_path, a_edge_to_ctg = load_tiling_paths(a_ctg_tiling_path, 'A');
    _, a_ctg_len = calc_tiling_paths_len(a_path)
    a_path = filter_tiling_paths_by_len(a_path, a_ctg_len, min_a_len)

    # Join the edge lookups.
    edge_to_ctg = {}
    edge_to_ctg.update(p_edge_to_ctg)
    edge_to_ctg.update(a_edge_to_ctg)

    # Join the tiling paths.
    all_tiling_paths = {}
    all_tiling_paths.update(p_path)
    all_tiling_paths.update(a_path)

    # Process the graph edges and extract only edges
    # and reads in the final graph.
    if tiling:
        read_in_graph, link_lines = get_gfa_links_from_tiling_paths(all_tiling_paths, edge_to_ctg);
    else:
        # Load the string graph.
        G_asm = AsmGraph(sg_edges_list, utg_data, ctg_paths)
        read_in_graph, link_lines = get_gfa_links_from_sg(G_asm, edge_to_ctg)

    # Load lengths of reads in the graph, and sequences if required.
    seq_bases, seq_len = load_seqs(preads_fasta, read_in_graph, write_reads)

    # Output version.
    fp_out.write(GFA_H_TAG + '\tVN:Z:1.0\n')

    # Output reads.
    for r in sorted(list(read_in_graph)):
        fp_out.write('\t'.join([GFA_S_TAG, r, GFA_SEQ_UNKNOWN if (not write_reads) else seq_bases[r], 'LN:i:%s' % seq_len[r]]) + '\n')

    # If required, output contigs.
    if (write_contigs):
        f = FastaReader(p_ctg_fasta)
        for r in f:
            if len(r.sequence) >= min_p_len:
                fp_out.write('\t'.join([GFA_S_TAG, r.name.split()[0], r.sequence, 'LN:i:%d' % (len(r.sequence))]) + '\n');
        f = FastaReader(a_ctg_fasta)
        for r in f:
            if len(r.sequence) >= min_a_len:
                fp_out.write('\t'.join([GFA_S_TAG, r.name.split()[0], r.sequence, 'LN:i:%d' % (len(r.sequence))]) + '\n');

    # Output links.
    for link_line in link_lines:
        fp_out.write(link_line + '\n')

    # Output contig paths.
    # Paths are already filtered by length, so just output them.
    for k in sorted(all_tiling_paths.keys()):
        out = format_gfa_path_line(k, all_tiling_paths[k], seq_len);
        fp_out.write('\t'.join(out) + '\n')

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Generates GFA output from FALCON's assembly.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--out', '-o', type=str, default='-', help='output GFA file to write to. If "-", stdout will be used')
    parser.add_argument('--p-ctg-tiling-path', type=str, default='p_ctg_tiling_path', help='location of the p_ctg tiling path file')
    parser.add_argument('--a-ctg-tiling-path', type=str, default='a_ctg_tiling_path', help='location of the a_ctg tiling path file')
    parser.add_argument('--preads-fasta', type=str, default='../1-preads_ovl/db2falcon/preads4falcon.fasta', help='path to the preads4falcon.fasta file')
    parser.add_argument('--p-ctg-fasta', type=str, default='p_ctg.fa', help='path to the primary contigs file')
    parser.add_argument('--a-ctg-fasta', type=str, default='a_ctg.fa', help='path to the associate contigs file')
    parser.add_argument('--sg-edges-list', type=str, default='sg_edges_list', help='string graph edges file from Falcon assembly')
    parser.add_argument('--utg-data', type=str, default='utg_data', help='unitig data file from Falcon')
    parser.add_argument('--ctg-paths', type=str, default='ctg_paths', help='contig paths file from Falcon assembly')
    parser.add_argument('--tiling', '-t', action='store_true', help="outputs only the tiling paths of contigs/associated contigs instead of the entire graph")
    parser.add_argument('--write-reads', '-r', action='store_true', help="output read sequences in S lines")
    parser.add_argument('--write-contigs', '-c', action='store_true', help="output contig sequences as S lines")
    parser.add_argument('--min-p-len', type=int, default=0, help='primary contig paths with length smaller than this will not be reported')
    parser.add_argument('--min-a-len', type=int, default=0, help='associate contig paths with length smaller than this will not be reported')
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)

    gfa_from_assembly(**vars(args))
