from falcon_kit.FastaReader import FastaReader
import argparse
import sys

def parse_args(argv):
    parser = argparse.ArgumentParser(description='remove duplicate a-tig, it assumes the working directory has the a_ctg_all.fa file')
    parser.add_argument('--max_idt', type=int, help="keep a-tig if the identity (in %) to the primary contig is <= max_idt", default = 96)
    parser.add_argument('--max_aln_cov', type=int, help="keep a-tig if the alignment coverage (in %) on the a-tig is <= max_aln_cov", default = 97)
    parser.add_argument('--min_len_diff', type=int, help="keep a-tig if the length different > min_len_diff", default = 500)
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    reads = FastaReader("a_ctg_all.fa")
    with open("a_ctg.fa","w") as f:
        for r in reads:
            tig_id, v, w, len_, ovl, ne, delta_l, idt, cov = r.name.split()
            if 100*float(idt) > args.max_idt and 100*float(cov) > args.max_aln_cov and\
               abs(int(delta_l)) < args.min_len_diff:
                   continue
            print >>f, ">"+r.name
            print >>f, r.sequence
