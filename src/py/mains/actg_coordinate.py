from falcon_kit.FastaReader import FastaReader


def main(argv=None):
  p_ctg_coor_map = {}
  with open("p_ctg_tiling_path") as f:
    for row in f:
        row = row.strip().split()
        ctg_id, v, w, edge_rid, b, e  = row[:6]
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


  a_ctg_fasta = FastaReader("a_ctg.fa")
  for r in a_ctg_fasta:
    rid = r.name.split()
    rid, v, w = rid[:3]
    pid = rid.split("-")[0]
    print rid, p_ctg_coor_map[pid][v], p_ctg_coor_map[pid][w]
