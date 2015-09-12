import networkx as nx
from falcon_kit.fc_asm_graph import AsmGraph


def main(argv=None):
  G_asm = AsmGraph("sg_edges_list", "utg_data", "ctg_paths")


  p_ctg_coor_map = {}
  for fn in ("p_ctg_tiling_path", "a_ctg_tiling_path"):
    f = open(fn)
    for row in f:
        row = row.strip().split()
        ctg_id, v, w, edge_rid, b, e  = row[:6]
        if ctg_id not in p_ctg_coor_map:
            coor = 0   # the p_ctg_tiling_path should be sorted by contig the order of the edges in the tiling path
            p_ctg_coor_map[ctg_id] = {}
            p_ctg_coor_map[ctg_id][v] = 0
            coor += abs(int(b) - int(e))
            p_ctg_coor_map[ctg_id][w] = coor
            G_asm.node_to_ctg[w]
            print ctg_id, v, 0, " ".join(list(G_asm.node_to_ctg[v]))
            print ctg_id, w, coor, " ".join(list(G_asm.node_to_ctg[w]))
            continue
        else:
            coor += abs(int(b) - int(e))
            p_ctg_coor_map[ctg_id][w] = coor 
            print ctg_id, w, coor, " ".join(list(G_asm.node_to_ctg[w]))
    f.close()
