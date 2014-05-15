from pbcore.io import FastaReader
import networkx as nx
import sys

u_graph = nx.DiGraph()
u_edges = {}
with open("./unit_edges.dat") as f:
    for l in f:
        v, w, path, seq = l.strip().split()
        u_edges.setdefault( (v, w), [] )
        u_edges[ (v, w) ].append( (path, seq) )
        u_graph.add_edge(v, w)
        
len(u_edges)
u_graph_r = u_graph.reverse()


p_tig_path = {}
a_tig_path = {}
with open("primary_tigs_paths_c") as f:
    for l in f:
        l = l.strip().split()
        id_ = l[0][1:]
        path = l[1:]
        p_tig_path[id_] = path

with open("all_tigs_paths") as f:
    for l in f:
        l = l.strip().split()
        id_ = l[0][1:]
        path = l[1:]
        a_tig_path[id_] = path

p_ugraph = nx.DiGraph()
p_sgraph = nx.DiGraph()
p_tig_id = sys.argv[1]

main_path = p_tig_path["%s_00" % p_tig_id]
all_nodes = set(main_path[:])
main_path_nodes = set(main_path[:])
p_ugraph.add_path(main_path)
for id_ in a_tig_path:
    if id_[:4] == p_tig_id:
        a_path = a_tig_path[id_]
        if a_path[0] in main_path_nodes and a_path[-1] in main_path_nodes:
            p_ugraph.add_path(a_path)
            for pp in a_path:
                all_nodes.add(pp)
        
for v, w in u_edges:
    if v in all_nodes and w in all_nodes:
        for p, s in u_edges[(v,w)]:
            p = p.split("-")
            p_sgraph.add_path(p)
            #print p
            for pp in p:
                all_nodes.add(pp)

nx.write_gexf(p_ugraph, "p_ugraph.gexf")
nx.write_gexf(p_sgraph, "p_sgraph.gexf")


preads = FastaReader(sys.argv[2])

all_nodes_ids = set( [s.split(":")[0] for s in list(all_nodes)] )
with open("p_sgraph_nodes.fa","w") as f:
    for r in preads:
        if r.name in all_nodes_ids:
            print >>f, ">"+r.name
            print >>f, r.sequence
