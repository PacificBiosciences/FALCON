import argparse
import os
import sys

from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader
from falcon_kit.gfa_graph import *

def run(fp_out, collected_gfa):
    with open(collected_gfa, 'r') as fp_in:
        gfa_graph = deserialize_gfa(fp_in)

    gfa_graph.write_gfa_v2(fp_out)

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Generates GFA output (on stdout) from FALCON's assembly.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('collected_gfa', type=str, default='asm.gfa.json',
                        help='Path to the file with collected and formatted data for generating the GFA')
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)

    run(sys.stdout, **vars(args))

if __name__ == '__main__':  # pragma: no cover
    main()
