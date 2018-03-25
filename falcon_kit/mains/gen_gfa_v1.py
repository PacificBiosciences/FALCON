from __future__ import absolute_import


from future.utils import viewitems
import argparse
import os
import sys

from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader
from falcon_kit.gfa_graph import GFAGraph


def calc_node_coords(tiling_path):
    """
    For a single tiling path (tiling_path is a list
    of edges for a particular contig) calculates the
    genomic coordinate of every node in the path.
    In case there are cycles in the tiling path,
    the existing node's coordinate will be overwritten.
    """
    if not tiling_path:
        return {}, 0
    coord_map = {}
    contig_len = 0
    v0, w0, b, e, l, idt, etype = tiling_path[0]
    coord_map[v0] = 0
    for v, w, b, e, l, idt, etype in tiling_path:
        if v not in coord_map:
            raise Exception(
                'Tiling path is not in sorted order. Node "{v!r}" does not yet have an assigned coordinate.'.format(v=v))
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
    for row in fp:
        row = row.strip().split()
        ctg_id, v, w, edge_rid, b, e, l, idt = row[:8]
        paths.setdefault(ctg_id, [])
        paths[ctg_id].append((v, w, int(b), int(e), int(l), float(idt), type_))
        edge_to_ctg[(v, w)] = ctg_id, type_
    return paths, edge_to_ctg


def load_tiling_paths(tiling_path_file, type_):
    """
    Open 'tiling_path_file' and call load_tiling_paths_from_stream().
    """
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
    for (ctg_id, edges) in viewitems(tiling_paths):
        path_coords[ctg_id], paths_len[ctg_id] = calc_node_coords(edges)
    return path_coords, paths_len


def filter_tiling_paths_by_len(tiling_paths, paths_len, min_len):
    """
    Given a dict of tiling_paths (key is contig name, value is a
    list of tiling path edges) creates and returns a new dict with
    tiling paths which are larger than the specified minimum length.
    """
    ret_paths = {}
    for (ctg_id, edges) in viewitems(tiling_paths):
        plen = paths_len[ctg_id]
        if plen >= min_len:
            ret_paths[ctg_id] = edges
    return ret_paths

def add_tiling_paths_to_gfa(p_ctg_fasta, a_ctg_fasta,
                            p_ctg_tiling_path, a_ctg_tiling_path,
                            min_p_len, min_a_len, gfa_graph, filter_tiling_paths_by_ctgid):
    # Associate tiling paths are not deduplicated.
    # We need the headers of the final haplotigs to filter
    # out the unnecessary tiling paths.
    a_ctg_headers = set()
    f = FastaReader(a_ctg_fasta)
    for r in f:
        a_ctg_headers.add(r.name.split(' ')[0])

    # Load and filter primary contig paths.
    p_paths, p_edge_to_ctg = load_tiling_paths(p_ctg_tiling_path, 'P')
    _, p_ctg_len = calc_tiling_paths_len(p_paths)
    p_paths = filter_tiling_paths_by_len(p_paths, p_ctg_len, min_p_len)
    p_paths = filter_tiling_paths_by_ctgid(p_paths)
    for (ctg_id, path) in viewitems(p_paths):
        gfa_graph.add_tiling_path(path, ctg_id)

    # Load and filter associate contig paths.
    a_paths, a_edge_to_ctg = load_tiling_paths(a_ctg_tiling_path, 'A')
    _, a_ctg_len = calc_tiling_paths_len(a_paths)
    a_paths = filter_tiling_paths_by_len(a_paths, a_ctg_len, min_a_len)
    a_paths = filter_tiling_paths_by_ctgid(a_paths)
    for (ctg_id, path) in viewitems(a_paths):
        if ctg_id in a_ctg_headers:
            gfa_graph.add_tiling_path(path, ctg_id)


def get_filter_tpbc(only_these_contigs):
    """Given a string, return either a no-op filter (when bool(str) is False)
    or a filter which keeps "only these contigs".
    """
    if only_these_contigs:
        # then engage an actual filter
        ctgs_to_include = set(open(only_these_contigs).read().splitlines())
        def filter_tiling_paths_by_ctgid(tiling_paths):
            """Filter out any contigs that don't exist in the set"""
            return {k:v for k,v in [x for x in viewitems(tiling_paths) if x[0].split('-')[0] in ctgs_to_include]}
        return filter_tiling_paths_by_ctgid
    def noop_filter(tiling_paths):
        return tiling_paths
    return noop_filter


def gfa_from_assembly(fp_out, p_ctg_tiling_path, a_ctg_tiling_path,
                      preads_fasta, p_ctg_fasta, a_ctg_fasta,
                      sg_edges_list, utg_data, ctg_paths,
                      add_string_graph, write_reads, write_contigs,
                      min_p_len, min_a_len, only_these_contigs):
    """
    This method produces the GFA-1 formatted output of the
    FALCON assembly.
    The graphical output is produced from either the entire string
    graph (only the non-filtered edges are considered) or from only
    the tiling paths. String graph can show the neighborhood of contig
    breaks, whereas the tiling path output is more sparse.
    Output is written to stdout.
    """
    gfa_graph = GFAGraph()

    filter_tiling_paths_by_ctgid = get_filter_tpbc(only_these_contigs)

    add_tiling_paths_to_gfa(p_ctg_fasta, a_ctg_fasta,
                            p_ctg_tiling_path, a_ctg_tiling_path,
                            min_p_len, min_a_len,
                            gfa_graph, filter_tiling_paths_by_ctgid)

    if add_string_graph:
        # Load the string graph.
        asm_graph = AsmGraph(sg_edges_list, utg_data, ctg_paths)
        gfa_graph.add_asm_graph(asm_graph)


    gfa_graph.write_gfa_v1(fp_out, preads_fasta, [
                           p_ctg_fasta, a_ctg_fasta], write_reads, write_contigs)


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
    parser.add_argument('--utg-data', type=str,
                        default='utg_data', help='unitig data file from Falcon')
    parser.add_argument('--ctg-paths', type=str, default='ctg_paths',
                        help='contig paths file from Falcon assembly')
    parser.add_argument('--add-string-graph', action='store_true',
                        help="in addition to tiling paths, output other edges and nodes from the final string graph")
    parser.add_argument('--write-reads', '-r', action='store_true',
                        help="output read sequences in S lines")
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

    gfa_from_assembly(sys.stdout, **vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
