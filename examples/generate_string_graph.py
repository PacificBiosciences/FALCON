from pbcore.io import FastaReader
import networkx as nx
import os
import shlex
import sys
import subprocess

class SGNode(object):
    def __init__(self, node_name):
        self.name = node_name
        self.out_edges = []
        self.in_edges = []
    def add_out_edge(self, out_edge):
        self.out_edges.append(out_edge)
    def add_in_edge(self, in_edge):
        self.in_edges.append(in_edge)

class SGEdge(object):
    def __init__(self, in_node, out_node):
        self.in_node = in_node
        self.out_node = out_node
        self.attr = {}
    def set_attribute(self, attr, value):
        self.attr[attr] = value

class StringGraph(object):
    def __init__(self):
        self.nodes = {}
        self.edges = {}
        self.n_mark = {}
        self.e_reduce = {}
        self.repeat_overlap = {}
        
    def add_node(self, node_name):
        if node_name not in self.nodes:
            self.nodes[node_name] = SGNode(node_name)
    
    def add_edge(self, in_node_name, out_node_name, **attributes):
        if (in_node_name, out_node_name) not in self.edges:
        
            self.add_node(in_node_name)
            self.add_node(out_node_name)
            in_node = self.nodes[in_node_name]
            out_node = self.nodes[out_node_name]    
            
            edge = SGEdge(in_node, out_node)
            self.edges[ (in_node_name, out_node_name) ] = edge
            in_node.add_out_edge(edge)
            out_node.add_in_edge(edge)
        edge =  self.edges[ (in_node_name, out_node_name) ]
        for k, v in attributes.items():
            edge.attr[k] = v
            
    def mark_tr_edges(self):
        n_mark = self.n_mark
        e_reduce = self.e_reduce
        FUZZ = 500
        for n in self.nodes:
            n_mark[n] = "vacant"
        for e in self.edges:
            e_reduce[e] = False
    
        for n_name, node in self.nodes.items():

            out_edges = node.out_edges
            if len(out_edges) == 0:
                continue
            
            out_edges.sort(key=lambda x: x.attr["length"])
            
            for e in out_edges:
                w = e.out_node
                n_mark[ w.name ] = "inplay"
            
            max_len = out_edges[-1].attr["length"]
            #longest_edge = out_edges[-1]
                
            max_len += FUZZ
            
            for e in out_edges:
                e_len = e.attr["length"]
                w = e.out_node
                if n_mark[w.name] == "inplay":
                    w.out_edges.sort( key=lambda x: x.attr["length"] )
                    for e2 in w.out_edges:
                        if e2.attr["length"] + e_len < max_len:
                            x = e2.out_node
                            if n_mark[x.name] == "inplay":
                                n_mark[x.name] = "eliminated"
            
            for e in out_edges:
                e_len = e.attr["length"]
                w = e.out_node
                w.out_edges.sort( key=lambda x: x.attr["length"] )
                if len(w.out_edges) > 0:
                    x = w.out_edges[0].out_node
                    if n_mark[x.name] == "inplay":
                        n_mark[x.name] = "eliminated"
                for e2 in w.out_edges:
                    if e2.attr["length"] < FUZZ:
                        x = e2.out_node
                        if n_mark[x.name] == "inplay":
                            n_mark[x.name] = "eliminated"
                            
            for out_edge in out_edges:
                v = out_edge.in_node
                w = out_edge.out_node
                if n_mark[w.name] == "eliminated":
                    e_reduce[ (v.name, w.name) ] = True
                n_mark[w.name] = "vacant"
                
    def mark_repeat_overlap(self):
        repeat_overlap = self.repeat_overlap
        in_degree = {}
        for n in self.nodes:
            c = 0
            for e in self.nodes[n].in_edges:
                v = e.in_node
                w = e.out_node
                if self.e_reduce[(v.name, w.name)] == False:
                    c += 1
            in_degree[n] = c
            #print n,c
        #print len([x for x in in_degree.items() if x[1]>1])
         
        for e_n, e in self.edges.items():
            v = e.in_node
            w = e.out_node
            if self.e_reduce[(v.name, w.name)] == False:
                repeat_overlap[ (v.name, w.name) ] = False
            else:
                repeat_overlap[ (v.name, w.name) ] = True
            
        for n in self.nodes:
            if len(self.nodes[n].out_edges) < 2:
                continue
            min_in_deg = None
            for e in self.nodes[n].out_edges:
                v = e.in_node
                w = e.out_node
                #print n, v.name, w.name
                if self.e_reduce[ (v.name, w.name) ] == True:
                    continue
                if min_in_deg == None:
                    min_in_deg = in_degree[w.name]
                    continue
                if in_degree[w.name] < min_in_deg:
                    min_in_deg = in_degree[w.name]
                #print n, w.name, in_degree[w.name]
            for e in self.nodes[n].out_edges:
                v = e.in_node
                w = e.out_node
                assert (v.name, w.name) in self.edges
                if in_degree[w.name] > min_in_deg:
                    if self.e_reduce[(v.name, w.name)] == False:
                        repeat_overlap[ (v.name, w.name) ] = True
                        
                    
        for e_n, e in self.edges.items():
            v = e.in_node
            w = e.out_node
            if repeat_overlap[ (v.name, w.name) ] == True:
                self.e_reduce[(v.name, w.name)] == True
                
                
    def get_out_edges_for_node(self, name, mask=True):
        rtn = []
        for e in self.nodes[name].out_edges:
            v = e.in_node
            w = e.out_node
            if self.e_reduce[ (v.name, w.name) ] == False:
                rtn.append(e)
        return rtn
        
        
    def get_in_edges_for_node(self, name, mask=True):
        rtn = []
        for e in self.nodes[name].in_edges:
            v = e.in_node
            w = e.out_node
            if self.e_reduce[ (v.name, w.name) ] == False:
                rtn.append(e)
        return rtn


RCMAP = dict(zip("ACGTacgtNn-","TGCAtgcaNn-"))
def generate_contig_from_path(sg, seqs, path):
    subseqs = []
    r_id, end = path[0].split(":")
    if end == "B":
        subseqs= [ "".join( [RCMAP[c] for c in seqs[r_id][::-1]] ) ]
    else:
        subseqs=[ seqs[r_id] ]
    
    count = 0
    for i in range( len( path ) -1 ):
        w_n, v_n = path[i:i+2]
        edge = sg.edges[ (w_n, v_n ) ]
        read_id, coor = edge.attr["label"].split(":")
        b,e = coor.split("-")
        b = int(b)
        e = int(e)
        if b < e:
            subseqs.append( seqs[read_id][b:e] )
        else:
            subseqs.append( "".join( [RCMAP[c] for c in seqs[read_id][b:e:-1]] ) )

    return "".join(subseqs)


def generate_contig(sg, seqs, out_fn=None):
    G = SGToNXG(sg)
    #out_fasta = open(out_fn, "w")
    count = 0
    out_f = open(out_fn, "w")
    out_f2 = open("edges.fa", "w")
    while len(G) > 0:
        root_nodes = set()
        for n in G:
            if G.in_degree(n) != 0 and G.out_degree(n) !=0 :
                root_nodes.add(n)
        if len(root_nodes) == 0: break
        candidates = []
        for n in list(root_nodes):
            sp =nx.single_source_shortest_path_length(G,n)
            sp = sp.items()
            sp.sort(key=lambda x : x[1])
            longest = sp[-1]
            candidates.append ( (longest[1], n, longest[0]) )
        candidates.sort()
        candidate = candidates[-1]
        path = nx.shortest_path(G, candidate[1], candidate[2])
        print >> out_f, ">%04d-%s-%s" % (count,path[0], path[-1])
        print >> out_f, generate_contig_from_path(sg, seqs, path)
        count += 1

        for i in range( len( path ) -1 ):
            w_n, v_n = path[i:i+2]
            edge = sg.edges[ (w_n, v_n ) ]
            read_id, coor = edge.attr["label"].split(":")
            b,e = coor.split("-")
            b = int(b)
            e = int(e)
            if b < e:
                print >> out_f2, ">"+ w_n + "-"+ v_n
                print >> out_f2, seqs[read_id][b:e] 
            else:
                print >> out_f2, ">"+ w_n + "-"+ v_n
                print >> out_f2, "".join( [RCMAP[c] for c in seqs[read_id][b:e:-1]] ) 

        n = path[-1]
        r_id, end = n.split(":")
        if end == "E":
            end = "B"
        elif end == "B":
            end = "E"
        if r_id+":"+end in G.nodes():
            compliment_nodes = list(nx.dfs_postorder_nodes(G, r_id+":"+end))
            for nn in compliment_nodes:
                G.remove_node(nn)
                      
        removed_count = 0
        
        for i in range( len( path ) -1 ):
            v_n, w_n = path[i:i+2]
            if ( v_n, w_n ) in G.edges():
                G.remove_edge( v_n, w_n )
        
        for n in path:
            if G.in_degree(n) == 0 and G.out_degree(n) == 0:
                G.remove_node(n)
                removed_count += 1            
        if removed_count == 0:
            break
    out_f.close()
    out_f2.close()
    nx.write_gml(G, "left_over.gml")





G=nx.Graph()
edges =set()
overlap_data = []
contained_reads = set()
with open("./out") as f:
    for l in f:
        l = l.strip().split()
        f_id, g_id, score, identity = l[:4]
        if f_id == g_id:
            continue
        score = int(score)
        identity = float(identity)
        contained = l[12]
        if contained == "contained":
            contained_reads.add(f_id)
            continue
        if contained == "contains":
            contained_reads.add(g_id)
            continue
        if contained == "none":
            continue
        if identity < 96:
            continue
        if score > -500:
            continue
        f_strain, f_start, f_end, f_len = (int(c) for c in l[4:8])
        g_strain, g_start, g_end, g_len = (int(c) for c in l[8:12])
        
        if f_start > 24 and f_len - f_end > 24:
            continue
        
        if g_start > 24 and g_len - g_end > 24:
            continue
        #if g_strain != 0:
        #    continue
        overlap_data.append( (f_id, g_id, score, identity,
                              f_strain, f_start, f_end, f_len,
                              g_strain, g_start, g_end, g_len) )


overlap_set = set()
sg = StringGraph()
#G=nx.Graph()
for od in overlap_data:
    f_id, g_id, score, identity = od[:4]
    if f_id in contained_reads:
        continue
    if g_id in contained_reads:
        continue
    f_s, f_b, f_e, f_l = od[4:8]
    g_s, g_b, g_e, g_l = od[8:12]
    overlap_pair = [f_id, g_id]
    overlap_pair.sort()
    overlap_pair = tuple( overlap_pair )
    if overlap_pair in overlap_set:
        continue
    else:
        overlap_set.add(overlap_pair)
    if g_s == 1:
        g_b, g_e = g_e, g_b
    if f_b > 10:
        if g_b < g_e:
            """
                 f.B         f.E
              f  ----------->
              g         ------------->
                        g.B           g.E
            """
            sg.add_edge( "%s:B" % g_id, "%s:B" % f_id, label = "%s:%d-%d" % (f_id, f_b, 0), 
                                                       length = abs(f_b-0),
                                                       score = -score)
            sg.add_edge( "%s:E" % f_id, "%s:E" % g_id, label = "%s:%d-%d" % (g_id, g_e, g_l), 
                                                       length = abs(g_e-g_l),
                                                       score = -score)
        else:
            """
                 f.B         f.E
              f  ----------->
              g         <-------------
                        g.E           g.B           
            """
            sg.add_edge( "%s:E" % g_id, "%s:B" % f_id, label = "%s:%d-%d" % (f_id, f_b, 0), 
                                                       length = abs(f_b -0),
                                                       score = -score)
            sg.add_edge( "%s:E" % f_id, "%s:B" % g_id, label = "%s:%d-%d" % (g_id, g_e, 0), 
                                                       length = abs(g_e- 0),
                                                       score = -score)
    else:
        if g_b < g_e:
            """
                                f.B         f.E
              f                 ----------->
              g         ------------->
                        g.B           g.E
            """
            sg.add_edge( "%s:B" % f_id, "%s:B" % g_id, label = "%s:%d-%d" % (g_id, g_b, 0), 
                                                       length = abs(g_b - 0),
                                                       score = -score)
            sg.add_edge( "%s:E" % g_id, "%s:E" % f_id, label = "%s:%d-%d" % (f_id, f_e, f_l), 
                                                       length = abs(f_e-f_l),
                                                       score = -score)
        else:
            """
                                f.B         f.E
              f                 ----------->
              g         <-------------
                        g.E           g.B           
            """
            sg.add_edge( "%s:B" % f_id, "%s:E" % g_id, label = "%s:%d-%d" % (g_id, g_b, g_l), 
                                                       length = abs(g_b - g_l),
                                                       score = -score)
            sg.add_edge( "%s:B" % g_id, "%s:E" % f_id, label = "%s:%d-%d" % (f_id, f_e, f_l), 
                                                       length = abs(f_e - f_l),
                                                       score = -score)
    
sg.mark_tr_edges()
print sum( [1 for c in sg.e_reduce.values() if c == True] )
print sum( [1 for c in sg.e_reduce.values() if c == False] )
sg.mark_repeat_overlap()
print sum( [1 for c in sg.repeat_overlap.values() if c == True] )
print sum( [1 for c in sg.repeat_overlap.values() if c == False] )
print len(sg.e_reduce), len(sg.repeat_overlap)

def SGToNXG(sg):
    G=nx.DiGraph()
    for v, w in sg.edges:
        if sg.e_reduce[(v, w)] != True:
        ##if 1:
            G.add_node( v )
            G.add_node( w )
            label = sg.edges[ (v, w) ].attr["label"]
            score = sg.edges[ (v, w) ].attr["score"]
            G.add_edge( v, w, label = label, weight = 0.001*score )
            #print in_node_name, out_node_name
    return G
nx.write_gml(SGToNXG(sg), "string_graph.gml")


seqs = {}
#f = FastaReader("pre_assembled_reads.fa")
f = FastaReader("sim_read.fa")
for r in f:
    seqs[r.name] = r.sequence.upper()


generate_contig(sg, seqs, out_fn="contigs.fa")
