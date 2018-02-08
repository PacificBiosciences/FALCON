from __future__ import print_function
from falcon_kit.FastaReader import open_fasta_reader
import argparse
import sys


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Removes duplicate a-tig, iff *all* conditions are violated. Assumes the working directory has the a_ctg_all.fa file, and produces a_ctg.fa',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--max_idt', type=int,
                        help="keep a-tig if the identity (in %) to the primary contig is <= max_idt", default=96)
    parser.add_argument('--max_aln_cov', type=int,
                        help="keep a-tig if the alignment coverage (in %) on the a-tig is <= max_aln_cov", default=97)
    parser.add_argument('--min_len_diff', type=int,
                        help="keep a-tig if the length different > min_len_diff", default=500)
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    with open_fasta_reader("a_ctg_all.fa") as reads:
        with open("a_ctg.fa", "w") as f:
            for r in reads:
                tig_id, v, w, len_, ovl, ne, delta_l, idt, cov = r.name.split()
                if 100 * float(idt) > args.max_idt and 100 * float(cov) > args.max_aln_cov and\
                   abs(int(delta_l)) < args.min_len_diff:
                    continue
                print(">" + r.name, file=f)
                print(r.sequence, file=f)


if __name__ == "__main__":
    main(sys.argv)
