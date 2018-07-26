"""
TODO: (from convo w/ Ivan)
the issue with this script (but would still like to re-read it to refresh my memory). The script loads all edge sequences and tries to do two things at once: create p_ctg and a_ctg sequences, and align the bubbles using those sequences


If we generate:
1. All paths first (as tiling paths) for all p_ctg and all a_ctg without loading sequences - this should not consume much space (take a look at *_tiling_paths files).
2. Load the first read of each tiling path fully, and only edge sequences for every transition, we can generate the output sequences with the same memory/disk consumption.
3. Align bubbles after that.

Our resource consumption should be same

Bubbles?
It aligns them to produce the identity score

After that the dedup_a_tigs.py script is used to deduplicate fake a_ctg.
But that script is simple, and only depends on the alignment info that the previous script stored in the a_ctg header.
"""
from __future__ import absolute_import
from __future__ import print_function


from builtins import zip
from builtins import range
import argparse
import sys
import networkx as nx
#from pbcore.io import FastaReader
from falcon_kit.FastaReader import open_fasta_reader
from falcon_kit import kup, falcon, DWA

read_fasta = "preads4falcon.fasta"
edge_data_file = "sg_edges_list"
utg_data_file = "utg_data"
ctg_data_file = "ctg_paths"

RCMAP = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))


class TooLongError(Exception): pass


def log(msg):
    sys.stderr.write(msg)
    sys.stderr.write('\n')


def rc(seq):
    return "".join([RCMAP[c] for c in seq[::-1]])


def get_aln_data(t_seq, q_seq):
    aln_data = []
    #x = []
    #y = []
    K = 8
    seq0 = t_seq
    lk_ptr = kup.allocate_kmer_lookup(1 << (K * 2))
    sa_ptr = kup.allocate_seq(len(seq0))
    sda_ptr = kup.allocate_seq_addr(len(seq0))
    kup.add_sequence(0, K, seq0, len(seq0), sda_ptr, sa_ptr, lk_ptr)
    q_id = "dummy"

    kmer_match_ptr = kup.find_kmer_pos_for_seq(
        q_seq, len(q_seq), K, sda_ptr, lk_ptr)
    kmer_match = kmer_match_ptr[0]

    if kmer_match.count != 0:
        aln_range_ptr = kup.find_best_aln_range(kmer_match_ptr, K, K * 5, 12)
        aln_range = aln_range_ptr[0]
        #x, y = list(zip(* [(kmer_match.query_pos[i], kmer_match.target_pos[i])
        #              for i in range(kmer_match.count)]))

        s1, e1, s2, e2 = aln_range.s1, aln_range.e1, aln_range.s2, aln_range.e2

        max_len = 250000 # to keep allocations < 16GB, given band_tol=1500
        if (e1 - s1) >= max_len or (e2 - s2) >= max_len:
            # DW.align() would crash, so raise here.
            # (500000 is the approx. upper bound for int overflow,
            #  but some users run out of memory anyway.)
            raise TooLongError('q_len={} or t_len={} are too big, over 500k'.format(
                (e1-s1), (e2-s2)))
        if e1 - s1 > 100:
            log('Calling DW_banded.align(q, {}, t, {}, 1500, 1)'.format(
                e1-s1, e2-s2))
            alignment = DWA.align(q_seq[s1:e1], e1 - s1,
                                  seq0[s2:e2], e2 - s2,
                                  1500, 1)

            if alignment[0].aln_str_size > 100:
                aln_data.append((q_id, 0, s1, e1, len(q_seq), s2, e2, len(
                    seq0), alignment[0].aln_str_size, alignment[0].dist))
                aln_str1 = alignment[0].q_aln_str
                aln_str0 = alignment[0].t_aln_str

            DWA.free_alignment(alignment)

        kup.free_aln_range(aln_range_ptr)

    kup.free_kmer_match(kmer_match_ptr)
    kup.free_kmer_lookup(lk_ptr)
    kup.free_seq_array(sa_ptr)
    kup.free_seq_addr_array(sda_ptr)
    return aln_data #, x, y


def reverse_end(node_id):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end


def yield_first_seq(one_path_edges, seqs):
            if one_path_edges and one_path_edges[0][0] != one_path_edges[-1][1]:
                # If non-empty, and non-circular,
                # prepend the entire first read.
                (vv, ww) = one_path_edges[0]
                (vv_rid, vv_letter) = vv.split(":")
                if vv_letter == 'E':
                    first_seq = seqs[vv_rid]
                else:
                    assert vv_letter == 'B'
                    first_seq = "".join([RCMAP[c] for c in seqs[vv_rid][::-1]])
                yield first_seq


def run(improper_p_ctg, proper_a_ctg):
    """improper==True => Neglect the initial read.
    We used to need that for unzip.
    """
    reads_in_layout = set()
    with open(edge_data_file) as f:
        for l in f:
            l = l.strip().split()
            """001039799:E 000333411:E 000333411 17524 20167 17524 99.62 G"""
            v, w, rid, s, t, aln_score, idt, type_ = l
            if type_ != "G":
                continue
            r1 = v.split(":")[0]
            reads_in_layout.add(r1)
            r2 = w.split(":")[0]
            reads_in_layout.add(r2)

    seqs = {}
    # load all p-read name into memory
    with open_fasta_reader(read_fasta) as f:
        for r in f:
            if r.name not in reads_in_layout:
                continue
            seqs[r.name] = r.sequence.upper() # name == rid-string

    edge_data = {}
    with open(edge_data_file) as f:
        for l in f:
            l = l.strip().split()
            """001039799:E 000333411:E 000333411 17524 20167 17524 99.62 G"""
            v, w, rid, s, t, aln_score, idt, type_ = l

            if type_ != "G":
                continue
            r1, dir1 = v.split(":")
            reads_in_layout.add(r1) # redundant, but harmless
            r2, dir2 = w.split(":")
            reads_in_layout.add(r2) # redundant, but harmless

            s = int(s)
            t = int(t)
            aln_score = int(aln_score)
            idt = float(idt)

            if s < t:
                e_seq = seqs[rid][s:t]
                assert 'E' == dir2
            else:
                # t and s were swapped for 'c' alignments in ovlp_to_graph.generate_string_graph():702
                # They were translated from reverse-dir to forward-dir coordinate system in LA4Falcon.
                e_seq = "".join([RCMAP[c] for c in seqs[rid][t:s][::-1]])
                assert 'B' == dir2
            edge_data[(v, w)] = (rid, s, t, aln_score, idt, e_seq)

    utg_data = {}
    with open(utg_data_file) as f:
        for l in f:
            l = l.strip().split()
            s, v, t, type_, length, score, path_or_edges = l
            if type_ not in ["compound", "simple", "contained"]:
                continue
            length = int(length)
            score = int(score)
            if type_ in ("simple", "contained"):
                path_or_edges = path_or_edges.split("~")
            else:
                path_or_edges = [tuple(e.split("~"))
                                 for e in path_or_edges.split("|")]
            utg_data[(s, v, t)] = type_, length, score, path_or_edges

    p_ctg_out = open("p_ctg.fa", "w")
    a_ctg_out = open("a_ctg_all.fa", "w")
    a_ctg_base_out = open("a_ctg_base.fa", "w")
    p_ctg_t_out = open("p_ctg_tiling_path", "w")
    a_ctg_t_out = open("a_ctg_tiling_path", "w")
    a_ctg_base_t_out = open("a_ctg_base_tiling_path", "w")
    layout_ctg = set()

    with open(ctg_data_file) as f:
        for l in f:
            l = l.strip().split()
            ctg_id, c_type_, i_utig, t0, length, score, utgs = l
            ctg_id = ctg_id
            s0 = i_utig.split("~")[0]

            if (reverse_end(t0), reverse_end(s0)) in layout_ctg:
                continue
            else:
                layout_ctg.add((s0, t0))

            ctg_label = i_utig + "~" + t0
            length = int(length)
            utgs = utgs.split("|")
            one_path = []
            total_score = 0
            total_length = 0

            #a_ctg_data = []
            a_ctg_group = {}

            for utg in utgs:
                s, v, t = utg.split("~")
                type_, length, score, path_or_edges = utg_data[(s, v, t)]
                total_score += score
                total_length += length
                if type_ == "simple":
                    if len(one_path) != 0:
                        one_path.extend(path_or_edges[1:])
                    else:
                        one_path.extend(path_or_edges)
                if type_ == "compound":

                    c_graph = nx.DiGraph()

                    all_alt_path = []
                    for ss, vv, tt in path_or_edges:
                        type_, length, score, sub_path = utg_data[(ss, vv, tt)]

                        v1 = sub_path[0]
                        for v2 in sub_path[1:]:
                            c_graph.add_edge(
                                v1, v2, e_score=edge_data[(v1, v2)][3])
                            v1 = v2

                    shortest_path = nx.shortest_path(c_graph, s, t, "e_score")
                    score = nx.shortest_path_length(c_graph, s, t, "e_score")
                    all_alt_path.append((score, shortest_path))

                    # a_ctg_data.append( (s, t, shortest_path) ) #first path is the same as the one used in the primary contig
                    while 1:
                        n0 = shortest_path[0]
                        for n1 in shortest_path[1:]:
                            c_graph.remove_edge(n0, n1)
                            n0 = n1
                        try:
                            shortest_path = nx.shortest_path(
                                c_graph, s, t, "e_score")
                            score = nx.shortest_path_length(
                                c_graph, s, t, "e_score")
                            #a_ctg_data.append( (s, t, shortest_path) )
                            all_alt_path.append((score, shortest_path))

                        except nx.exception.NetworkXNoPath:
                            break
                        # if len(shortest_path) < 2:
                        #    break
                    all_alt_path.sort()
                    all_alt_path.reverse()
                    shortest_path = all_alt_path[0][1]
                    if len(one_path) != 0:
                        one_path.extend(shortest_path[1:])
                    else:
                        one_path.extend(shortest_path)

                    a_ctg_group[(s, t)] = all_alt_path

            if len(one_path) == 0:
                continue

            one_path_edges = list(zip(one_path[:-1], one_path[1:]))

            if improper_p_ctg:
                sub_seqs = []
            else:
                sub_seqs = list(yield_first_seq(one_path_edges, seqs))
            for vv, ww in one_path_edges:
                rid, s, t, aln_score, idt, e_seq = edge_data[(vv, ww)]
                sub_seqs.append(e_seq)
                print("%s %s %s %s %d %d %d %0.2f" % (
                    ctg_id, vv, ww, rid, s, t, aln_score, idt), file=p_ctg_t_out)
            print(">%s %s %s %d %d" % (
                ctg_id, ctg_label, c_type_, total_length, total_score), file=p_ctg_out)
            print("".join(sub_seqs), file=p_ctg_out)

            a_id = 1
            for v, w, in a_ctg_group:
                # get the base sequence used in the primary contig
                atig_output = []

                score, atig_path = a_ctg_group[(v, w)][0]
                atig_path_edges = list(zip(atig_path[:-1], atig_path[1:]))
                if not proper_a_ctg:
                    sub_seqs = []
                else:
                    sub_seqs = list(yield_first_seq(atig_path_edges, seqs))
                total_length = 0
                total_score = 0
                for vv, ww in atig_path_edges:
                    rid, s, t, aln_score, idt, e_seq = edge_data[(vv, ww)]
                    sub_seqs.append(e_seq)
                    total_length += abs(s - t)
                    total_score += aln_score

                base_seq = "".join(sub_seqs)
                atig_output.append(
                    (v, w, atig_path, total_length, total_score, base_seq, atig_path_edges, 0, 1, 1))

                for score, atig_path in a_ctg_group[(v, w)][1:]:
                    atig_path_edges = list(zip(atig_path[:-1], atig_path[1:]))
                    if not proper_a_ctg:
                        sub_seqs = []
                    else:
                        sub_seqs = list(yield_first_seq(atig_path_edges, seqs))
                    total_length = 0
                    total_score = 0
                    for vv, ww in atig_path_edges:
                        rid, s, t, aln_score, idt, e_seq = edge_data[(vv, ww)]
                        sub_seqs.append(e_seq)
                        total_length += abs(s - t)
                        total_score += aln_score

                    seq = "".join(sub_seqs)

                    delta_len = len(seq) - len(base_seq)
                    idt = 0.0
                    cov = 0.0
                    if len(base_seq) > 2000 and len(seq) > 2000:
                        try:
                          aln_data = get_aln_data(base_seq, seq)
                          if len(aln_data) != 0:
                            idt = 1.0 - 1.0 * \
                                aln_data[-1][-1] / aln_data[-1][-2]
                            cov = 1.0 * \
                                (aln_data[-1][3] - aln_data[-1]
                                 [2]) / aln_data[-1][4]
                        except TooLongError:
                            log('WARNING: Seqs were too long for get_aln_data(), so we set idt/cov low enough to prevent filtering by dedup_a_tigs, at atig_path[:-1] == {}'.format(atig_path[:-1]))
                            idt = -1.0
                            cov = -1.0

                    atig_output.append(
                        (v, w, atig_path, total_length, total_score, seq, atig_path_edges, delta_len, idt, cov))

                if len(atig_output) == 1:
                    continue

                sub_id = 0
                for data in atig_output:
                    v0, w0, tig_path, total_length, total_score, seq, atig_path_edges, delta_len, a_idt, cov = data
                    for vv, ww in atig_path_edges:
                        rid, s, t, aln_score, idt, e_seq = edge_data[(vv, ww)]
                        if sub_id != 0:
                            print("%s-%03d-%02d %s %s %s %d %d %d %0.2f" % (
                                ctg_id, a_id, sub_id, vv, ww, rid, s, t, aln_score, idt), file=a_ctg_t_out)
                        else:
                            print("%s-%03d-%02d %s %s %s %d %d %d %0.2f" % (
                                ctg_id, a_id, sub_id, vv, ww, rid, s, t, aln_score, idt), file=a_ctg_base_t_out)

                    if sub_id != 0:
                        print(">%s-%03d-%02d %s %s %d %d %d %d %0.2f %0.2f" % (
                            ctg_id, a_id, sub_id, v0, w0, total_length, total_score, len(atig_path_edges), delta_len, a_idt, cov), file=a_ctg_out)
                        print(seq, file=a_ctg_out)
                    else:
                        print(">%s-%03d-%02d %s %s %d %d %d %d %0.2f %0.2f" % (
                            ctg_id, a_id, sub_id, v0, w0, total_length, total_score, len(atig_path_edges), delta_len, a_idt, cov), file=a_ctg_base_out)
                        print(seq, file=a_ctg_base_out)

                    sub_id += 1

                a_id += 1

    a_ctg_out.close()
    a_ctg_base_out.close()
    p_ctg_out.close()
    a_ctg_t_out.close()
    a_ctg_base_t_out.close()
    a_ctg_t_out.close()
    p_ctg_t_out.close()


def main(argv=sys.argv):
    description = 'Generate the primary and alternate contig fasta files and tiling paths, given the string graph.'
    epilog = """
We assume these input files, in cwd:

    read_fasta = "preads4falcon.fasta"
    edge_data_file = "sg_edges_list"
    utg_data_file = "utg_data"
    ctg_data_file = "ctg_paths"

We write these:

    p_ctg_out = open("p_ctg.fa", "w")
    a_ctg_out = open("a_ctg_all.fa", "w")
    a_ctg_base_out = open("a_ctg_base.fa", "w")
    p_ctg_t_out = open("p_ctg_tiling_path", "w")
    a_ctg_t_out = open("a_ctg_tiling_path", "w")
    a_ctg_base_t_out = open("a_ctg_base_tiling_path", "w")
"""
    parser = argparse.ArgumentParser(
            description=description,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=epilog)
    parser.add_argument('--improper-p-ctg', action='store_true',
            help='Skip the initial read in each p_ctg path.')
    parser.add_argument('--proper-a-ctg', action='store_true',
            help='Skip the initial read in each a_ctg path.')
    args = parser.parse_args(argv[1:])
    run(**vars(args))


if __name__ == "__main__":
    main(sys.argv)
