from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

from builtins import range
from ctypes import (POINTER, c_char_p, c_uint, c_uint,
                    c_uint, c_uint, c_uint, c_double, string_at)
from falcon_kit.multiproc import Pool
from falcon_kit import falcon
import argparse
import os
import re
import sys
import falcon_kit


falcon.generate_consensus.argtypes = [
    POINTER(c_char_p), c_uint, c_uint, c_uint, c_double]
falcon.generate_consensus.restype = POINTER(falcon_kit.ConsensusData)
falcon.free_consensus_data.argtypes = [POINTER(falcon_kit.ConsensusData)]


def get_longest_reads(seqs, max_n_read, max_cov_aln, sort=True):
    # including the sort kwarg allows us to avoid a redundant sort
    # in get_consensus_trimmed()
    if sort:
        seqs = seqs[:1] + sorted(seqs[1:], key=lambda x: -len(x))

    longest_n_reads = max_n_read
    if max_cov_aln > 0:
        longest_n_reads = 1
        seed_len = len(seqs[0])
        read_cov = 0
        for seq in seqs[1:]:
            if read_cov // seed_len > max_cov_aln:
                break
            longest_n_reads += 1
            read_cov += len(seq)

        longest_n_reads = min(longest_n_reads, max_n_read)

    return(seqs[:longest_n_reads])


def get_alignment(seq1, seq0, edge_tolerance=1000):

    kup = falcon_kit.kup
    K = 8
    lk_ptr = kup.allocate_kmer_lookup(1 << (K * 2))
    sa_ptr = kup.allocate_seq(len(seq0))
    sda_ptr = kup.allocate_seq_addr(len(seq0))
    kup.add_sequence(0, K, seq0, len(seq0), sda_ptr, sa_ptr, lk_ptr)

    kup.mask_k_mer(1 << (K * 2), lk_ptr, 16)
    kmer_match_ptr = kup.find_kmer_pos_for_seq(
        seq1, len(seq1), K, sda_ptr, lk_ptr)
    kmer_match = kmer_match_ptr[0]
    aln_range_ptr = kup.find_best_aln_range2(kmer_match_ptr, K, K * 50, 25)
    #x,y = zip( * [ (kmer_match.query_pos[i], kmer_match.target_pos[i]) for i in range(kmer_match.count )] )
    aln_range = aln_range_ptr[0]
    kup.free_kmer_match(kmer_match_ptr)
    s1, e1, s0, e0, km_score = aln_range.s1, aln_range.e1, aln_range.s2, aln_range.e2, aln_range.score
    e1 += K + K // 2
    e0 += K + K // 2
    kup.free_aln_range(aln_range)
    len_1 = len(seq1)
    len_0 = len(seq0)
    if e1 > len_1:
        e1 = len_1
    if e0 > len_0:
        e0 = len_0

    aln_size = 1
    if e1 - s1 > 500:

        aln_size = max(e1 - s1, e0 - s0)
        aln_score = int(km_score * 48)
        aln_q_s = s1
        aln_q_e = e1
        aln_t_s = s0
        aln_t_e = e0

    kup.free_seq_addr_array(sda_ptr)
    kup.free_seq_array(sa_ptr)
    kup.free_kmer_lookup(lk_ptr)

    if s1 > edge_tolerance and s0 > edge_tolerance:
        return 0, 0, 0, 0, 0, 0, "none"

    if len_1 - e1 > edge_tolerance and len_0 - e0 > edge_tolerance:
        return 0, 0, 0, 0, 0, 0, "none"

    if e1 - s1 > 500 and aln_size > 500:
        return s1, s1 + aln_q_e - aln_q_s, s0, s0 + aln_t_e - aln_t_s, aln_size, aln_score, "aln"
    else:
        return 0, 0, 0, 0, 0, 0, "none"


def get_consensus_without_trim(c_input):
    seqs, seed_id, config = c_input
    min_cov, K, max_n_read, min_idt, edge_tolerance, trim_size, min_cov_aln, max_cov_aln = config
    if len(seqs) > max_n_read:
        seqs = get_longest_reads(seqs, max_n_read, max_cov_aln, sort=True)
    seqs_ptr = (c_char_p * len(seqs))()
    seqs_ptr[:] = seqs
    consensus_data_ptr = falcon.generate_consensus(
        seqs_ptr, len(seqs), min_cov, K, min_idt)

    consensus = string_at(consensus_data_ptr[0].sequence)[:]
    eff_cov = consensus_data_ptr[0].eff_cov[:len(consensus)]
    falcon.free_consensus_data(consensus_data_ptr)
    del seqs_ptr
    return consensus, seed_id


def get_consensus_with_trim(c_input):
    seqs, seed_id, config = c_input
    min_cov, K, max_n_read, min_idt, edge_tolerance, trim_size, min_cov_aln, max_cov_aln = config
    trim_seqs = []
    seed = seqs[0]
    for seq in seqs[1:]:
        aln_data = get_alignment(seq, seed, edge_tolerance)
        s1, e1, s2, e2, aln_size, aln_score, c_status = aln_data
        if c_status == "none":
            continue
        if aln_score > 1000 and e1 - s1 > 500:
            e1 -= trim_size
            s1 += trim_size
            trim_seqs.append((e1 - s1, seq[s1:e1]))
    trim_seqs.sort(key=lambda x: -x[0])  # use longest alignment first
    trim_seqs = [x[1] for x in trim_seqs]

    trim_seqs = [seed] + trim_seqs
    if len(trim_seqs[1:]) > max_n_read:
        # seqs already sorted, dont' sort again
        trim_seqs = get_longest_reads(
            trim_seqs, max_n_read, max_cov_aln, sort=False)

    seqs_ptr = (c_char_p * len(trim_seqs))()
    seqs_ptr[:] = trim_seqs
    consensus_data_ptr = falcon.generate_consensus(
        seqs_ptr, len(trim_seqs), min_cov, K, min_idt)
    consensus = string_at(consensus_data_ptr[0].sequence)[:]
    eff_cov = consensus_data_ptr[0].eff_cov[:len(consensus)]
    falcon.free_consensus_data(consensus_data_ptr)
    del seqs_ptr
    return consensus, seed_id


def get_seq_data(config, min_n_read, min_len_aln):
    max_len = 100000
    min_cov, K, max_n_read, min_idt, edge_tolerance, trim_size, min_cov_aln, max_cov_aln = config
    seqs = []
    seed_id = None
    seed_len = 0
    seqs_data = []
    read_cov = 0
    read_ids = set()
    with sys.stdin as f:
        for l in f:
            l = l.strip().split()
            if len(l) != 2:
                continue

            read_id = l[0]
            seq = l[1]
            if len(seq) > max_len:
                seq = seq[:max_len - 1]

            if read_id not in ("+", "-", "*"):
                if len(seq) >= min_len_aln:
                    if len(seqs) == 0:
                        seqs.append(seq)  # the "seed"
                        seed_len = len(seq)
                        seed_id = read_id
                    if read_id not in read_ids:  # avoidng using the same read twice. seed is used again here by design
                        seqs.append(seq)
                        read_ids.add(read_id)
                        read_cov += len(seq)
            elif l[0] == "+":
                if len(seqs) >= min_n_read and read_cov // seed_len >= min_cov_aln:
                    seqs = get_longest_reads(
                        seqs, max_n_read, max_cov_aln, sort=True)
                    yield (seqs, seed_id, config)
                #seqs_data.append( (seqs, seed_id) )
                seqs = []
                read_ids = set()
                seed_id = None
                read_cov = 0
            elif l[0] == "*":
                seqs = []
                read_ids = set()
                seed_id = None
                read_cov = 0
            elif l[0] == "-":
                # yield (seqs, seed_id)
                #seqs_data.append( (seqs, seed_id) )
                break


def format_seq(seq, col):
    return "\n".join([seq[i:(i + col)] for i in range(0, len(seq), col)])


def main(argv=sys.argv):
    parser = argparse.ArgumentParser(description='a simple multi-processor consensus sequence generator',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--n_core', type=int, default=24,
                        help='number of processes used for generating consensus; '
                        '0 for main process only')
    parser.add_argument('--min_cov', type=int, default=6,
                        help='minimum coverage to break the consensus')
    parser.add_argument('--min_cov_aln', type=int, default=10,
                        help='minimum coverage of alignment data; a seed read with less than MIN_COV_ALN average depth' +
                        ' of coverage will be completely ignored')
    parser.add_argument('--max_cov_aln', type=int, default=0,  # 0 to emulate previous behavior
                        help='maximum coverage of alignment data; a seed read with more than MAX_COV_ALN average depth' + \
                        ' of coverage of the longest alignments will be capped, excess shorter alignments will be ignored')
    parser.add_argument('--min_len_aln', type=int, default=0,  # 0 to emulate previous behavior
                        help='minimum length of a sequence in an alignment to be used in consensus; any shorter sequence will be completely ignored')
    parser.add_argument('--min_n_read', type=int, default=10,
                        help='1 + minimum number of reads used in generating the consensus; a seed read with fewer alignments will ' +
                        'be completely ignored')
    parser.add_argument('--max_n_read', type=int, default=500,
                        help='1 + maximum number of reads used in generating the consensus')
    parser.add_argument('--trim', action="store_true", default=False,
                        help='trim the input sequence with k-mer spare dynamic programming to find the mapped range')
    parser.add_argument('--output_full', action="store_true", default=False,
                        help='output uncorrected regions too')
    parser.add_argument('--output_multi', action="store_true", default=False,
                        help='output multi correct regions')
    parser.add_argument('--min_idt', type=float, default=0.70,
                        help='minimum identity of the alignments used for correction')
    parser.add_argument('--edge_tolerance', type=int, default=1000,
                        help='for trimming, the there is unaligned edge leng > edge_tolerance, ignore the read')
    parser.add_argument('--trim_size', type=int, default=50,
                        help='the size for triming both ends from initial sparse aligned region')
    good_region = re.compile("[ACGT]+")
    args = parser.parse_args(argv[1:])

    def Start():
        print('Started a worker in %d from parent %d' % (
            os.getpid(), os.getppid()), file=sys.stderr)
    exe_pool = Pool(args.n_core, initializer=Start)
    if args.trim:
        get_consensus = get_consensus_with_trim
    else:
        get_consensus = get_consensus_without_trim

    K = 8
    config = args.min_cov, K, \
        args.max_n_read, args.min_idt, args.edge_tolerance, args.trim_size, args.min_cov_aln, args.max_cov_aln
    # TODO: pass config object, not tuple, so we can add fields
    for res in exe_pool.imap(get_consensus, get_seq_data(config, args.min_n_read, args.min_len_aln)):
        cns, seed_id = res
        if len(cns) < 500:
            continue

        if args.output_full:
            print(">" + seed_id + "_f")
            print(cns)
        else:
            cns = good_region.findall(cns)
            if len(cns) == 0:
                continue
            if args.output_multi:
                seq_i = 0
                for cns_seq in cns:
                    if len(cns_seq) < 500:
                        continue
                    if seq_i >= 10:
                        break
                    print(">prolog/%s%01d/%d_%d" % (seed_id, seq_i, 0, len(cns_seq)))
                    print(format_seq(cns_seq, 80))
                    seq_i += 1
            else:
                cns.sort(key=lambda x: len(x))
                print(">" + seed_id)
                print(cns[-1])


if __name__ == "__main__":
    main(sys.argv)
