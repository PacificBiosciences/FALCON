#!/usr/bin/env python

#################################################################################$$
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$

from pbcore.io import FastaReader
import networkx as nx
import os
import shlex
import sys
import subprocess

DEBUG_LOG_LEVEL = 0

class SGNode(object):
    """
    class representing a node in the string graph
    """
    def __init__(self, node_name):
        self.name = node_name
        self.out_edges = []
        self.in_edges = []
    def add_out_edge(self, out_edge):
        self.out_edges.append(out_edge)
    def add_in_edge(self, in_edge):
        self.in_edges.append(in_edge)

class SGEdge(object):
    """
    class representing an edge in the string graph
    """
    def __init__(self, in_node, out_node):
        self.in_node = in_node
        self.out_node = out_node
        self.attr = {}
    def set_attribute(self, attr, value):
        self.attr[attr] = value

def reverse_end( node_id ):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end

class StringGraph(object):
    """
    class representing the string graph
    """
    def __init__(self):
        self.nodes = {}
        self.edges = {}
        self.n_mark = {}
        self.e_reduce = {}
        self.repeat_overlap = {}
        
    def add_node(self, node_name):
        """ 
        add a node into the graph by given a node name
        """
        if node_name not in self.nodes:
            self.nodes[node_name] = SGNode(node_name)
    
    def add_edge(self, in_node_name, out_node_name, **attributes):
        """ 
        add an edge into the graph by given a pair of nodes
        """
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

    def init_reduce_dict(self):
        for e in self.edges:
            self.e_reduce[e] = False

    def mark_chimer_edge(self):

        for e_n, e in self.edges.items():
            v = e_n[0]
            w = e_n[1]
            overlap_count = 0
            for w_out_e in self.nodes[w].out_edges:
                w_out_n = w_out_e.out_node.name
                if (v, w_out_n) in self.edges:
                    overlap_count += 1
            for v_in_e in self.nodes[v].in_edges:
                v_in_n = v_in_e.in_node.name
                if (v_in_n, w) in self.edges:
                    overlap_count += 1
            if self.e_reduce[ (v, w) ] != True:
                if overlap_count == 0:
                    self.e_reduce[(v, w)] = True
                    #print "XXX: chimer edge %s %s removed" % (v, w)
                    v, w = reverse_end(w), reverse_end(v)
                    self.e_reduce[(v, w)] = True
                    #print "XXX: chimer edge %s %s removed" % (v, w)



    def mark_spur_edge(self):

        for  v in self.nodes:
            if len(self.nodes[v].out_edges) > 1:
                for out_edge in self.nodes[v].out_edges:
                    w = out_edge.out_node.name
                    
                    if len(self.nodes[w].out_edges) == 0 and self.e_reduce[ (v, w) ] != True:
                        #print "XXX: spur edge %s %s removed" % (v, w)
                        self.e_reduce[(v, w)] = True
                        v2, w2 = reverse_end(w), reverse_end(v)
                        #print "XXX: spur edge %s %s removed" % (v2, w2)
                        self.e_reduce[(v, w)] = True

            if len(self.nodes[v].in_edges) > 1:
                for in_edge in self.nodes[v].in_edges:
                    w = in_edge.in_node.name
                    if len(self.nodes[w].in_edges) == 0 and self.e_reduce[ (w, v) ] != True:
                        #print "XXX: spur edge %s %s removed" % (w, v)
                        self.e_reduce[(w, v)] = True
                        v2, w2 = reverse_end(w), reverse_end(v)
                        #print "XXX: spur edge %s %s removed" % (w2, v2)
                        self.e_reduce[(w, v)] = True


    def mark_tr_edges(self):
        """
        transitive reduction
        """
        n_mark = self.n_mark
        e_reduce = self.e_reduce
        FUZZ = 500
        for n in self.nodes:
            n_mark[n] = "vacant"
    
        for n_name, node in self.nodes.items():

            out_edges = node.out_edges
            if len(out_edges) == 0:
                continue
            
            out_edges.sort(key=lambda x: x.attr["length"])
            
            for e in out_edges:
                w = e.out_node
                n_mark[ w.name ] = "inplay"
            
            max_len = out_edges[-1].attr["length"]
                
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
                    #print "XXX: tr edge %s %s removed" % (v.name, w.name)
                    v_name, w_name = reverse_end(w.name), reverse_end(v.name)
                    e_reduce[(v_name, w_name)] = True
                    #print "XXX: tr edge %s %s removed" % (v_name, w_name)
                n_mark[w.name] = "vacant"
                

    def mark_best_overlap(self):
        """
        find the best overlapped edges
        """

        best_edges = set()

        for v in self.nodes:

            out_edges = self.nodes[v].out_edges
            if len(out_edges) > 0:
                out_edges.sort(key=lambda e: e.attr["score"])
                e = out_edges[-1]
                best_edges.add( (e.in_node.name, e.out_node.name) )

            in_edges = self.nodes[v].in_edges
            if len(in_edges) > 0:
                in_edges.sort(key=lambda e: e.attr["score"])
                e = in_edges[-1]
                best_edges.add( (e.in_node.name, e.out_node.name) )

        if DEBUG_LOG_LEVEL > 1:
            print "X", len(best_edges)

        for e_n, e in self.edges.items():
            v = e_n[0]
            w = e_n[1]
            if self.e_reduce[ (v, w) ] != True:
                if (v, w) not in best_edges:
                    self.e_reduce[(v, w)] = True
                    #print "XXX: in best edge %s %s removed" % (v, w)
                    v2, w2 = reverse_end(w), reverse_end(v)
                    #print "XXX: in best edge %s %s removed" % (v2, w2)
                    self.e_reduce[(v2, w2)] = True
                
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

    def get_best_out_edge_for_node(self, name, mask=True):
        rtn = []
        for e in self.nodes[name].out_edges:
            v = e.in_node
            w = e.out_node
            if self.e_reduce[ (v.name, w.name) ] == False:
                rtn.append(e)
        rtn.sort(key=lambda e: e.attr["score"])

        return rtn[-1]

    def get_best_in_edge_for_node(self, name, mask=True):
        rtn = []
        for e in self.nodes[name].in_edges:
            v = e.in_node
            w = e.out_node
            if self.e_reduce[ (v.name, w.name) ] == False:
                rtn.append(e)
        rtn.sort(key=lambda e: e.attr["score"])
        return rtn[-1]
        

RCMAP = dict(zip("ACGTacgtNn-","TGCAtgcaNn-"))
def generate_seq_from_path(sg, seqs, path):
    subseqs = []
    r_id, end = path[0].split(":")
    
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


def reverse_path( path ):
    new_path = []
    for n in list(path[::-1]):
        rid, end = n.split(":")
        new_end = "B" if end == "E" else "E"
        new_path.append( rid+":"+new_end)
    return new_path


def generate_unitig(sg, seqs, out_fn, connected_nodes = None):

    """
    given a string graph:sg and the sequences: seqs, write the unitig fasta file into out_fn
    the funtion return a reduct graph representing the reduce string graph where the edges are unitigs
    
    some extra files generated: 
        unit_edges.dat : an easy to parse file for unitig data
        unit_edge_paths : the file contains the information of the path of all unitigs
        uni_graph.gexf: the unitig graph in gexf format for visulization
    """

    G = SGToNXG(sg)
    if connected_nodes != None:
        connected_nodes = set(sg.nodes)
    out_fasta = open(out_fn, "w")
    nodes_for_tig = set()
    sg_edges = set()
    for v, w in sg.edges:
        if sg.e_reduce[(v, w)] != True:
            sg_edges.add( (v, w) )
    count = 0
    edges_in_tigs = set()

    uni_edges = {}
    path_f = open("unit_edge_paths","w")
    uni_edge_f = open("unit_edges.dat", "w")
    while len(sg_edges) > 0:
        v, w = sg_edges.pop()

        #nodes_for_tig.remove(n)
        upstream_nodes = []
        
        c_node = v
        p_in_edges = sg.get_in_edges_for_node(c_node)
        p_out_edges = sg.get_out_edges_for_node(c_node)
        while len(p_in_edges) == 1 and len(p_out_edges) == 1:
            p_node = p_in_edges[0].in_node
            upstream_nodes.append(p_node.name)
            if (p_node.name, c_node) not in  sg_edges:
                break
            p_in_edges = sg.get_in_edges_for_node(p_node.name)
            p_out_edges = sg.get_out_edges_for_node(p_node.name)
            c_node = p_node.name

        upstream_nodes.reverse()  
            
        downstream_nodes = []
        c_node = w 
        n_out_edges = sg.get_out_edges_for_node(c_node)
        n_in_edges = sg.get_in_edges_for_node(c_node)
        while len(n_out_edges) == 1 and len(n_in_edges) == 1:
            n_node = n_out_edges[0].out_node
            downstream_nodes.append(n_node.name)
            if (c_node, n_node.name) not in  sg_edges:
                break
            n_out_edges = sg.get_out_edges_for_node(n_node.name)
            n_in_edges = sg.get_in_edges_for_node(n_node.name)
            c_node = n_node.name 
        
        whole_path = upstream_nodes + [v, w] + downstream_nodes
        count += 1
        subseq = generate_seq_from_path(sg, seqs, whole_path) 
        #subseq = ""
        uni_edges.setdefault( (whole_path[0], whole_path[-1]), [] )
        uni_edges[(whole_path[0], whole_path[-1])].append(  ( whole_path, subseq ) )
        print >> uni_edge_f, whole_path[0], whole_path[-1], "-".join(whole_path), subseq
        print >>path_f, ">%05dc-%s-%s-%d %s" % (count, whole_path[0], whole_path[-1], len(whole_path), " ".join(whole_path))
        print >>out_fasta, ">%05dc-%s-%s-%d" % (count, whole_path[0], whole_path[-1], len(whole_path))
        print >>out_fasta, subseq
        for i in range( len( whole_path ) -1 ):
            w_n, v_n = whole_path[i:i+2]
            try:
                sg_edges.remove( (w_n, v_n) )
            except KeyError: #if an edge is already deleted, ignore it
                pass

        r_whole_path = reverse_path( whole_path )
        count += 1
        subseq = generate_seq_from_path(sg, seqs, r_whole_path) 
        #subseq = ""
        uni_edges.setdefault( (r_whole_path[0], r_whole_path[-1]), [] )
        uni_edges[(r_whole_path[0], r_whole_path[-1])].append(  ( r_whole_path, subseq ) )
        print >> uni_edge_f, r_whole_path[0], r_whole_path[-1], "-".join(r_whole_path), subseq
        print >>path_f, ">%05dc-%s-%s-%d %s" % (count, r_whole_path[0], r_whole_path[-1], len(r_whole_path), " ".join(r_whole_path))
        print >>out_fasta, ">%05dc-%s-%s-%d" % (count, r_whole_path[0], r_whole_path[-1], len(r_whole_path))
        print >>out_fasta, subseq
        for i in range( len( r_whole_path ) -1 ):
            w_n, v_n = r_whole_path[i:i+2]
            try:
                sg_edges.remove( (w_n, v_n) )
            except KeyError: #if an edge is already deleted, ignore it
                pass


    path_f.close()
    uni_edge_f.close()
    #uni_graph = nx.DiGraph()
    #for n1, n2 in uni_edges.keys():
    #    uni_graph.add_edge(n1, n2, count = len( uni_edges[ (n1,n2) ] ))
    #nx.write_gexf(uni_graph, "uni_graph.gexf")

    out_fasta.close()
    return uni_edges

def neighbor_bound(G, v, w, radius):
    """
    test if the node v and the node w are connected within a radius in graph G
    """
    g1 = nx.ego_graph(G, v, radius=radius, undirected=False)
    g2 = nx.ego_graph(G, w, radius=radius, undirected=False)
    if len(set(g1.edges()) & set(g2.edges())) > 0:
        return True
    else:
        return False


def is_branch_node(G, n):
    """
    test whether the node n is a "branch node" which the paths from any of two of 
    its offsprings do not intersect within a given radius
    """
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


def get_bundle( path, u_graph, u_graph_r ):

    """ 
    find a sub-graph contain the nodes between the start and the end of the path
    inputs: 
        u_graph : a unitig graph
    returns:
        bundle_graph: the whole bundle graph 
        bundle_paths: the paths in the bundle graph 
        sub_graph2_edges: all edges of the bundle graph
    
    """

    p_start, p_end = path[0], path[-1]
    p_nodes = set(path)
    p_edges = set(zip(path[:-1], path[1:]))

    down_path = nx.ego_graph(u_graph, p_start, radius=len(p_nodes), undirected=False)
    up_path = nx.ego_graph(u_graph_r, p_end, radius=len(p_nodes), undirected=False)
    subgraph_nodes = set(down_path) & set(up_path)
    

    sub_graph = nx.DiGraph()
    for v, w in u_graph.edges_iter():
        if v in subgraph_nodes and w in subgraph_nodes:            
            if (v, w) in p_edges:
                sub_graph.add_edge(v, w, color = "red")
            else:
                sub_graph.add_edge(v, w, color = "black")

    sub_graph2 = nx.DiGraph()
    tips = set()
    tips.add(path[0])
    sub_graph_r = sub_graph.reverse()
    visited = set()
    ct = 0
    is_branch = is_branch_node(sub_graph, path[0]) #if the start node is a branch node
    if is_branch:
        n = tips.pop()
        e = sub_graph.out_edges([n])[0] #pick one path the build the subgraph
        sub_graph2.add_edge(e[0], e[1], n_weight = u_graph[e[0]][e[1]]["n_weight"])
        if e[1] not in visited:
            last_node = e[1]
            visited.add(e[1])
            r_id, orientation = e[1].split(":")
            orientation = "E" if orientation == "B" else "E"
            visited.add( r_id +":" + orientation)
            if not is_branch_node(sub_graph_r, e[1]): 
                tips.add(e[1])
        
    while len(tips) != 0:
        n = tips.pop()
        out_edges = sub_graph.out_edges([n])
        if len(out_edges) == 1:
            e = out_edges[0]
            sub_graph2.add_edge(e[0], e[1], n_weight = u_graph[e[0]][e[1]]["n_weight"])
            last_node = e[1]
            if e[1] not in visited:                       
                visited.add(e[1])
                r_id, orientation = e[1].split(":")
                orientation = "E" if orientation == "B" else "E"
                visited.add( r_id +":" + orientation)
                if not is_branch_node(sub_graph_r, e[1]): 
                    tips.add(e[1])
        else:
        
            is_branch = is_branch_node(sub_graph, n)
            if not is_branch:
                for e in out_edges:
                    sub_graph2.add_edge(e[0], e[1], n_weight = u_graph[e[0]][e[1]]["n_weight"])
                    last_node = e[1]
                    if e[1] not in visited:
                        r_id, orientation = e[1].split(":")
                        visited.add(e[1])
                        orientation = "E" if orientation == "B" else "E"
                        visited.add( r_id +":" + orientation)
                        if not is_branch_node(sub_graph_r, e[1]):
                            tips.add(e[1])
        ct += 1
    last_node = None
    longest_len = 0
        
    sub_graph2_nodes = sub_graph2.nodes()
    sub_graph2_edges = sub_graph2.edges()


    new_path = [path[0]]
    for n in sub_graph2_nodes:
        if len(sub_graph2.out_edges(n)) == 0 :
            path_t = nx.shortest_path(sub_graph2, source = path[0], target = n, weight = "n_weight")
            path_len = len(path_t)
            if path_len > longest_len:
                last_node = n
                longest_len = path_len
                new_path = path_t

    if last_node == None:
        for n in sub_graph2_nodes:
            path_t = nx.shortest_path(sub_graph2, source = path[0], target = n, weight = "n_weight")
            path_len = len(path_t)
            if path_len > longest_len:
                last_node = n
                longest_len = path_len
                new_path = path_t


    path = new_path

    # clean up sub_graph2 according to new begin and end
    sub_graph2_r = sub_graph2.reverse()
    down_path = nx.ego_graph(sub_graph2, path[0], radius=len(path), undirected=False)
    up_path = nx.ego_graph(sub_graph2_r, path[-1], radius=len(path), undirected=False)
    subgraph_nodes = set(down_path) & set(up_path)
    for v in sub_graph2_nodes:
        if v not in subgraph_nodes:
            sub_graph2.remove_node(v)
    
    if DEBUG_LOG_LEVEL > 1:
        print "new_path", path[0], last_node, len(sub_graph2_nodes), path


    bundle_paths = [path]
    p_nodes = set(path)
    p_edges = set(zip(path[:-1], path[1:]))

    sub_graph2_nodes = sub_graph2.nodes()
    sub_graph2_edges = sub_graph2.edges()

    nodes_idx = dict( [ (n[1], n[0]) for n in enumerate(path) ]  )
    
         
    # create a list of subpath that has no branch
    non_branch_subpaths = []
    wi = 0
    vi = 0
    v = path[0]
    while v != path[-1] and wi < len(path)-1:
        wi += 1
        w = path[wi]
        while len( sub_graph2.successors(w) ) == 1 and len( sub_graph2.predecessors(w) ) == 1 and wi < len(path)-1:
            wi += 1
            w = path[wi]
        if  len( sub_graph2.successors(v) )!= 1 or len( sub_graph2.predecessors(w) )!= 1:
            branched = True
        else:
            branched = False

        if not branched:
            non_branch_subpaths.append( path[vi:wi+1] )
        v = w
        vi = wi

    # create the accompany_graph that has the path of the alternative subpaths
    
    associate_graph = nx.DiGraph()
    for v, w in sub_graph2.edges_iter():
        if (v, w) not in p_edges:
            associate_graph.add_edge(v, w, n_weight = sub_graph2[v][w]["n_weight"])

    if DEBUG_LOG_LEVEL > 1:
        print "associate_graph size:", len(associate_graph)           
        print "non_branch_subpaths",len(non_branch_subpaths), non_branch_subpaths

    # construct the bundle graph                
    associate_graph_nodes = set(associate_graph.nodes())
    bundle_graph = nx.DiGraph()
    bundle_graph.add_path( path )
    for i in range(len(non_branch_subpaths)-1):
        if len(non_branch_subpaths[i]) == 0 or len( non_branch_subpaths[i+1] ) == 0:
            continue
        e1, e2 = non_branch_subpaths[i: i+2]
        v = e1[-1]
        w = e2[0]
        if v == w:
            continue
        in_between_node_count = nodes_idx[w] - nodes_idx[v] 
        if v in associate_graph_nodes and w in associate_graph_nodes:
            try:
                a_path = nx.shortest_path(associate_graph, v, w, "n_weight")    
            except nx.NetworkXNoPath:
                continue
            bundle_graph.add_path( a_path )      
            bundle_paths.append( a_path )

    return bundle_graph, bundle_paths, sub_graph2_edges
            
def get_bundles(u_edges):
    
    """
    input: all unitig edges
    output: the assembled primary_tigs.fa and all_tigs.fa
    """

    ASM_graph = nx.DiGraph()
    out_f = open("primary_tigs.fa", "w")
    main_tig_paths = open("primary_tigs_paths","w")
    sv_tigs = open("all_tigs.fa","w")
    sv_tig_paths = open("all_tigs_paths","w")
    max_weight = 0 
    for v, w in u_edges:
        x = max( [len(s[1]) for s in u_edges[ (v,w) ] ] )
        if DEBUG_LOG_LEVEL > 1:
            print "W", v, w, x
        if x > max_weight:
            max_weight = x
            
    in_edges = {}
    out_edges = {}
    for v, w in u_edges:
        in_edges.setdefault(w, []) 
        out_edges.setdefault(w, []) 
        in_edges[w].append( (v, w) )

        out_edges.setdefault(v, [])
        in_edges.setdefault(v, [])
        out_edges[v].append( (v, w) )

    u_graph = nx.DiGraph()
    for v,w in u_edges:

        u_graph.add_edge(v, w, n_weight = max_weight - max( [len(s[1]) for s in  u_edges[ (v,w) ] ] ) )
    
    bundle_edge_out = open("bundle_edges","w")
    bundle_index = 0


    components = nx.weakly_connected_component_subgraphs(u_graph)
    components = [ (len(c), c) for c in components ]
    components.sort()
    #components.reverse()
    allS = len(u_graph)
    ssG = 0.0
    processed_overlaps = set()
    for sG, G in components:

        ssG += sG
        print "process graph of size ", sG, "%0.2f %0.2f" % (ssG, ssG/allS)
        G_edges = set(G.edges())

        dual_component = False
        
        for v, w in list(G_edges):
            v = v.split(":")[0]
            w = w.split(":")[0]
            if (v, w) in processed_overlaps:
                dual_component = True
                break

        if dual_component == True:
            continue

        for v, w in list(G_edges):
            v = v.split(":")[0]
            w = w.split(":")[0]
            processed_overlaps.add( (v,w) )
            processed_overlaps.add( (w,v) )

        G_r = G.reverse()
        visited_u_edges = set()

        while len(G) > 0:
            out_f.flush()
            main_tig_paths.flush()
            sv_tigs.flush()
            sv_tig_paths.flush()
            
            
            #root_nodes = set() 
            candidates = [] 
            for n in G: 
                sp =nx.single_source_shortest_path_length(G, n) 
                sp = sp.items() 
                sp.sort(key=lambda x : x[1]) 
                longest = sp[-1] 
                if DEBUG_LOG_LEVEL > 2:
                    print "L", n, longest[0]
                if longest[0].split(":")[0] == n.split(":")[0]: #avoid a big loop 
                    continue
                candidates.append ( (longest[1], n, longest[0]) ) 

                n = longest[0]
                sp =nx.single_source_shortest_path_length(G_r, n) 
                sp = sp.items() 
                sp.sort(key=lambda x : x[1]) 
                longest = sp[-1] 
                if DEBUG_LOG_LEVEL > 2:
                    print "L", n, longest[0]
                if longest[0].split(":")[0] == n.split(":")[0]: #avoid a big loop 
                    continue
                candidates.append ( (longest[1], longest[0], n) ) 
                if len(candidates) != 0:
                    break

            if len(candidates) == 0:
                print "no more candiate", len(G.edges()), len(G.nodes())
                if len(G_edges) > 0:
                    path = G_edges.pop()
                    G_edges.add(path)
                    print path
                else:
                    break
            else:
                candidates.sort() 
                
                candidate = candidates[-1] 
                
                if candidate[1] == candidate[2]: 
                    G.remove_node(candidate[1]) 
                    G_r.remove_node(candidate[1])
                    continue 
             
                path = nx.shortest_path(G, candidate[1], candidate[2], "n_weight") 

            if DEBUG_LOG_LEVEL > 1:
                print "X", path[0], path[-1], len(path)
            
            cmp_edges = set()
            #g_edges = set(G.edges())
            new_path = []  
            tail = True
            # avioid confusion due to long palindrome sequence
            if len(path) > 2:
                for i in range( 0, len( path ) - 1 ):
                    v_n, w_n = path[i:i+2]
                    new_path.append(v_n)
                    # the comment out code below might be useful for filter out some high connectivity nodes
                    #if (v_n, w_n) in cmp_edges or\
                    #    len(u_graph.out_edges(w_n)) > 5 or\
                    #    len(u_graph.in_edges(w_n)) > 5:
                    if (v_n, w_n) in cmp_edges: 
                        tail = False
                        break

                    r_id, end = v_n.split(":")
                    end = "E" if end == "B" else "B" 
                    v_n2 = r_id + ":" + end 

                    r_id, end = w_n.split(":")
                    end = "E" if end == "B" else "B" 
                    w_n2 = r_id + ":" + end 

                    if (w_n2, v_n2) in G_edges:
                        cmp_edges.add( (w_n2, v_n2) )

                if tail:
                    new_path.append(w_n)
            else:
                new_path = path[:]
                    
            
            if len(new_path) > 1:
                path = new_path
                
                if DEBUG_LOG_LEVEL > 2:
                    print "Y", path[0], path[-1], len(path)

                bundle_graph, bundle_paths, bundle_graph_edges = get_bundle( path, G, G_r )
                for bg_edge in bundle_graph_edges:
                    print >> bundle_edge_out, bundle_index, "edge", bg_edge[0], bg_edge[1]
                for path_ in bundle_paths:
                    print >>bundle_edge_out, "path", bundle_index, " ".join(path_) 

                edges_to_be_removed = set()
                if DEBUG_LOG_LEVEL > 2:
                    print "Z", bundle_paths[0][0], bundle_paths[0][-1]
                    print bundle_index, len(path), len(bundle_paths[0]), len(bundle_paths), len(bundle_graph_edges)

                if len(bundle_graph_edges) > 0:

                    ASM_graph.add_path(bundle_paths[0], ctg="%04d" % bundle_index)
                    extra_u_edges = []
                    
                    print >> main_tig_paths, ">%04d %s" % ( bundle_index, " ".join(bundle_paths[0]) )
                    subseqs = []
                
                    for i in range(len(bundle_paths[0]) - 1): 
                        v, w = bundle_paths[0][i:i+2]
                        edges_to_be_removed.add( (v,w) )
                        uedges = u_edges[ (v,w) ]
                        uedges.sort( key= lambda x: len(x[0]) )
                        subseqs.append( uedges[-1][1] )
                        visited_u_edges.add( "-".join(uedges[-1][0]) ) 
                        for ue in uedges:
                            if "-".join(ue[0]) not in visited_u_edges:
                                visited_u_edges.add("-".join(ue[0]))
                                extra_u_edges.append(ue)
                    seq = "".join(subseqs)        
                    sv_tig_idx = 0
                    print >> sv_tig_paths, ">%04d-%04d %s" % ( bundle_index, sv_tig_idx, " ".join(bundle_paths[0]) )
                    if len(seq) > 0:
                        print >> out_f, ">%04d %s-%s" % (bundle_index, bundle_paths[0][0], bundle_paths[0][-1])
                        print >> out_f, seq
                        print >> sv_tigs, ">%04d-%04d %s-%s" % (bundle_index, sv_tig_idx, bundle_paths[0][0], bundle_paths[0][-1])
                        print >> sv_tigs, "".join(subseqs)

                    sv_tig_idx += 1

                    for sv_path in bundle_paths[1:]:
                        print >> sv_tig_paths, ">%04d-%04d %s" % ( bundle_index, sv_tig_idx, " ".join(sv_path) )
                        ASM_graph.add_path(sv_path, ctg="%04d" % bundle_index)
                        subseqs = []
                        for i in range(len(sv_path) - 1): 
                            v, w = sv_path[i:i+2]
                            edges_to_be_removed.add( (v,w) )
                            uedges = u_edges[ (v,w) ]
                            uedges.sort( key= lambda x: len(x[0]) )
                            subseqs.append( uedges[-1][1] )
                            visited_u_edges.add( "-".join(uedges[-1][0]) ) 
                            for ue in uedges:
                                if "-".join(ue[0]) not in visited_u_edges:
                                    visited_u_edges.add("-".join(ue[0]))
                                    extra_u_edges.append(ue)
                        seq = "".join(subseqs)        
                        if len(seq) > 0: 
                            print >> sv_tigs, ">%04d-%04d %s-%s" % (bundle_index, sv_tig_idx, sv_path[0], sv_path[-1])
                            print >> sv_tigs, "".join(subseqs)
                        sv_tig_idx += 1
                    for u_path, seq in extra_u_edges:
                        #u_path = u_path.split("-")
                        ASM_graph.add_edge(u_path[0], u_path[-1], ctg="%04d" % bundle_index)
                        print >> sv_tig_paths, ">%04d-%04d-u %s" % ( bundle_index, sv_tig_idx, " ".join(u_path) )
                        print >> sv_tigs, ">%04d-%04d-u %s-%s" % (bundle_index, sv_tig_idx, u_path[0], u_path[-1])
                        print >> sv_tigs, seq
                        sv_tig_idx += 1
                        
                    
                    bundle_index += 1
            else:
                #TODO, consolidate code here
                v,w = path
                uedges = u_edges[ (v,w) ]
                uedges.sort( key= lambda x: len(x[0]) )
                subseqs.append( uedges[-1][1] )
                seq = "".join(subseqs)
                print >> sv_tig_paths, ">%04d-%04d %s" % ( bundle_index, sv_tig_idx, " ".join(paths) )
                print >> sv_tigs, ">%04d-%04d-u %s-%s" % (bundle_index, sv_tig_idx, path[0], path[-1])
                print >> sv_tigs, seq
                sv_tig_idx += 1
                bundle_index += 1
                bundle_graph_edges = zip(path[:-1],path[1:])
            
            #clean up the graph

            edges = set(G.edges())
            edges_to_be_removed |= set(bundle_graph_edges)

            if DEBUG_LOG_LEVEL > 2:
                print "BGE",bundle_graph_edges
            
            edge_remove_count = 0
            for v, w in edges_to_be_removed:
                if (v, w) in edges:
                    G.remove_edge( v, w )
                    G_r.remove_edge( w, v )
                    G_edges.remove( (v, w) )
                    edge_remove_count += 1
                    if DEBUG_LOG_LEVEL > 2:
                        print "remove edge", bundle_index, w, v
                    
            edges = set(G.edges())
            for v, w in edges_to_be_removed:

                r_id, end = v.split(":")
                end = "E" if end == "B" else "B"
                v = r_id + ":" + end

                r_id, end = w.split(":")
                end = "E" if end == "B" else "B"
                w = r_id + ":" + end

                if (w, v) in edges:
                    G.remove_edge( w, v )
                    G_edges.remove( (w, v) )
                    G_r.remove_edge( v, w )
                    edge_remove_count += 1
                    if DEBUG_LOG_LEVEL > 2:
                        print "remove edge", bundle_index, w, v

            if edge_remove_count == 0:
                break
                
            nodes = G.nodes()
            for n in nodes:
                if G.in_degree(n) == 0 and G.out_degree(n) == 0:
                    G.remove_node(n)
                    G_r.remove_node(n)
                    if DEBUG_LOG_LEVEL > 2:
                        print "remove node", n 

    sv_tig_paths.close()
    sv_tigs.close()
    main_tig_paths.close()
    out_f.close()
    bundle_edge_out.close()
    return ASM_graph



def SGToNXG(sg):
    G=nx.DiGraph()

    max_score = max([ sg.edges[ e ].attr["score"] for e in sg.edges if sg.e_reduce[e] != True ])
    out_f = open("edges_list","w")
    for v, w in sg.edges:
        if sg.e_reduce[(v, w)] != True:
        ##if 1:
            out_degree = len(sg.nodes[v].out_edges)
            G.add_node( v, size = out_degree )
            G.add_node( w, size = out_degree )
            label = sg.edges[ (v, w) ].attr["label"]
            score = sg.edges[ (v, w) ].attr["score"]
            print >>out_f, v, w, label, score 
            G.add_edge( v, w, label = label, weight = 0.001*score, n_weight = max_score - score )
            #print in_node_name, out_node_name
    out_f.close()
    return G

if __name__ == "__main__":

    import argparse
    
    parser = argparse.ArgumentParser(description='a example string graph assembler that is desinged for handling diploid genomes')
    parser.add_argument('overlap_file', help='a file that contains the overlap information.')
    parser.add_argument('read_fasta', help='the file that contains the sequence to be assembled')
    parser.add_argument('--min_len', type=int, default=4000, 
                        help='minimum length of the reads to be considered for assembling')
    parser.add_argument('--min_idt', type=float, default=96,
                        help='minimum alignment identity of the reads to be considered for assembling')
    parser.add_argument('--disable_chimer_prediction', action="store_true", default=False,
                        help='you may want to disable this as some reads can be falsely identified as chimers in low coverage case')

    args = parser.parse_args()


    overlap_file = args.overlap_file
    read_fasta = args.read_fasta

    contained_reads = set()
    chimer_ids = set()

    with open("rc_out_all") as f:
        for l in f:
            l = l.strip().split()
            if l[1] == "2":
                chimer_ids.add(l[0])
            if l[1] == "1":
                contained_reads.add(l[0])
    print len(chimer_ids)
    
    seqs = {}
    # load all p-reads into memory
    f = FastaReader(read_fasta)
    for r in f:
        if r.name in contained_reads:
            continue
        if r.name in chimer_ids:
            continue
        seqs[r.name] = r.sequence.upper()

    G=nx.Graph()
    edges =set()
    overlap_data = []
    contained_reads = set()
    overlap_count = {}


    # loop through the overlapping data to load the data in the a python array
    # contained reads are identified 

    with open(overlap_file) as f:
        for l in f:
            l = l.strip().split()

            #work around for some ill formed data recored
            if len(l) != 13:
                continue
            
            f_id, g_id, score, identity = l[:4]
            if f_id == g_id:  # don't need self-self overlapping
                continue

            if g_id not in seqs: 
                continue

            if f_id not in seqs:
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

            if identity < args.min_idt: # only take record with >96% identity as overlapped reads
                continue
            #if score > -2000:
            #    continue
            f_strain, f_start, f_end, f_len = (int(c) for c in l[4:8])
            g_strain, g_start, g_end, g_len = (int(c) for c in l[8:12])

            # only used reads longer than the 4kb for assembly
            if f_len < args.min_len: continue
            if g_len < args.min_len: continue
            
            # double check for proper overlap
            if f_start > 24 and f_len - f_end > 24:  # allow 24 base tolerance on both sides of the overlapping
                continue
            
            if g_start > 24 and g_len - g_end > 24:
                continue
            
            if g_strain == 0:
                if f_start < 24 and g_len - g_end > 24:
                    continue
                if g_start < 24 and f_len - f_end > 24:
                    continue
            else:
                if f_start < 24 and g_start > 24:
                    continue
                if g_start < 24 and f_start > 24:
                    continue

            overlap_data.append( (f_id, g_id, score, identity,
                                  f_strain, f_start, f_end, f_len,
                                  g_strain, g_start, g_end, g_len) )

            overlap_count[f_id] = overlap_count.get(f_id,0)+1
            overlap_count[g_id] = overlap_count.get(g_id,0)+1
            
    print "###", len(overlap_data), len(contained_reads)
    overlap_set = set()
    sg = StringGraph()
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
        if overlap_pair in overlap_set:  # don't allow duplicated records
            continue
        else:
            overlap_set.add(overlap_pair)

        
        if g_s == 1: # revered alignment, swapping the begin and end coordinates
            g_b, g_e = g_e, g_b
        
        # build the string graph edges for each overlap
        if f_b > 24:
            if g_b < g_e:
                """
                     f.B         f.E
                  f  ----------->
                  g         ------------->
                            g.B           g.E
                """
                if f_b == 0 or g_e - g_l == 0:
                    continue
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
                if f_b == 0 or g_e == 0:
                    continue
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
                if g_b == 0 or f_e - f_l == 0:
                    continue
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
                if g_b - g_l == 0 or f_e - f_l ==0:
                    continue
                sg.add_edge( "%s:B" % f_id, "%s:E" % g_id, label = "%s:%d-%d" % (g_id, g_b, g_l), 
                                                           length = abs(g_b - g_l),
                                                           score = -score)
                sg.add_edge( "%s:B" % g_id, "%s:E" % f_id, label = "%s:%d-%d" % (f_id, f_e, f_l), 
                                                           length = abs(f_e - f_l),
                                                           score = -score)


    sg.init_reduce_dict()
    #if not args.disable_chimer_prediction:
    #    sg.mark_chimer_edge()
    sg.mark_spur_edge()
    sg.mark_tr_edges() # mark those edges that transitive redundant

    #if DEBUG_LOG_LEVEL > 1:
    if 1:
        print sum( [1 for c in sg.e_reduce.values() if c == True] )
        print sum( [1 for c in sg.e_reduce.values() if c == False] )

    sg.mark_best_overlap() # mark those edges that are best overlap edges

    if DEBUG_LOG_LEVEL > 1:
        print sum( [1 for c in sg.e_reduce.values() if c == False] )


    G = SGToNXG(sg)
    nx.write_gexf(G, "string_graph.gexf") # output the raw string string graph for visuliation
    nx.write_adjlist(G, "string_graph.adj") # write out the whole adjacent list of the string graph

    u_edges = generate_unitig(sg, seqs, out_fn = "unitgs.fa") # reduct to string graph into unitig graph
    ASM_graph = get_bundles(u_edges )  # get the assembly
    nx.write_gexf(ASM_graph, "asm_graph.gexf")
