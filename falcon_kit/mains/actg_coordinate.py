from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
from falcon_kit.FastaReader import open_fasta_reader
import sys


def main(argv=sys.argv):
    p_ctg_coor_map = {}
    with open("p_ctg_tiling_path") as f:
        for row in f:
            row = row.strip().split()
            ctg_id, v, w, edge_rid, b, e = row[:6]
            if ctg_id not in p_ctg_coor_map:
                coor = 0   # the p_ctg_tiling_path should be sorted by contig the order of the edges in the tiling path
                p_ctg_coor_map[ctg_id] = {}
                p_ctg_coor_map[ctg_id][v] = 0
                coor += abs(int(b) - int(e))
                p_ctg_coor_map[ctg_id][w] = coor
                continue
            else:
                coor += abs(int(b) - int(e))
                p_ctg_coor_map[ctg_id][w] = coor

    with open_fasta_reader("a_ctg.fa") as a_ctg_fasta:
        for r in a_ctg_fasta:
            rid = r.name.split()
            rid, v, w = rid[:3]
            pid = rid.split("-")[0]
            print(rid, p_ctg_coor_map[pid][v], p_ctg_coor_map[pid][w])


if __name__ == "__main__":
    main(sys.argv)
