import networkx as nx
from pbcore.io import FastaReader

def neighbor_bound(G, v, w, radius):
    g1 = nx.ego_graph(G, v, radius=radius, undirected=False)
    g2 = nx.ego_graph(G, w, radius=radius, undirected=False)
    if len(set(g1.edges()) & set(g2.edges())) > 0:
        return True
    else:
        return False
    
def is_branch_node(G, n):
    out_edges = G.out_edges([n])
    n2 = [ e[1] for e in out_edges ]
    is_branch = False
    for i in range(len(n2)):
        for j in range(i+1, len(n2)):
            v = n2[i]
            w = n2[j]
            if neighbor_bound(G, v, w, 10) == False:
                is_branch = True
                break
        if is_branch == True:
            break
    return is_branch


def get_seq(u_edges, path):
    subseqs = []
    for i in range( len(path) - 1):
        v, w = path[i:i+2]
        uedges = u_edges[ (v, w) ]
        uedges.sort( key= lambda x: len(x[0]) )
        subseqs.append( uedges[-1][1] )
        
    return "".join(subseqs)




u_edges = {}
with open("unit_edges.dat") as f:
    for l in f:
        v, w, path, seq = l.strip().split()
        u_edges.setdefault( (v, w), [] )
        u_edges[ (v, w) ].append( (path, seq) )
len(u_edges)


primary_tigs_path = {}
primary_path_graph = nx.DiGraph()
begin_nodes = {}
end_nodes ={}
with open("primary_tigs_paths") as f:
    for l in f:
        l = l.strip().split()
        name = l[0][1:]
        path = l[1:]
        primary_tigs_path[name] = path
        if len(path) < 3:
            continue
        for i in range(len(path)-1):
            n1 = path[i].split(":")[0]
            n2 = path[i+1].split(":")[0]
            primary_path_graph.add_edge( n1, n2)
        begin_nodes.setdefault(path[0], [])
        begin_nodes[path[0]].append( name )
        end_nodes.setdefault(path[-1], [])
        end_nodes[path[-1]].append( name )



path_names = primary_tigs_path.keys()
path_names.sort()
primary_path_graph_r = primary_path_graph.reverse()
path_f = open("primary_tigs_paths_c","w")
with open("primary_tigs_c.fa","w") as out_f:
    for name in path_names:
        sub_idx = 0
        c_path = [ primary_tigs_path[name][0] ]
        for v in primary_tigs_path[name][1:]:
            break_path = False
            
            vn = v.split(":")[0]

            if primary_path_graph.out_degree(vn) > 1:
                break_path = is_branch_node(primary_path_graph, vn)
            if primary_path_graph.in_degree(vn) > 1:
                break_path = is_branch_node(primary_path_graph_r, vn)
            if break_path:
                c_path.append(v)
                seq = get_seq(u_edges, c_path)
                print >>out_f, ">%s_%02d" % (name, sub_idx)
                print >>out_f, seq
                print >>path_f, ">%s_%02d" % (name, sub_idx), " ".join(c_path)
                #print c_path
                c_path = [v]
                sub_idx += 1
            else:
                c_path.append(v)
                
        if len(c_path) > 1:
            seq = get_seq(u_edges, c_path)
            print >>out_f, ">%s_%02d" % (name, sub_idx)
            print >>out_f, seq
            print >>path_f, ">%s_%02d" % (name, sub_idx), " ".join(c_path)

            
path_f.close()
