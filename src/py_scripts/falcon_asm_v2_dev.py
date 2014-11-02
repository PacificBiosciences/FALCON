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

#from pbcore.io import FastaReader
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

    def mark_chimer_edges(self):

        for n_name in self.nodes:
            n = self.nodes[n_name]
            
            out_nodes = set( [ e.out_node for e in n.out_edges ] )
            in_nodes = [e.in_node for e in n.in_edges ] 
            is_chimer = True
            for in_node in in_nodes:
                for v in [e.out_node for e in in_node.out_edges]:
                    if v in out_nodes:
                        is_chimer = False
                        break

            if is_chimer == True:
                for e in n.out_edges:
                    v, w =  e.in_node.name, e.out_node.name
                    self.e_reduce[ (v, w) ] = True
                for e in n.in_edges:
                    v, w =  e.in_node.name, e.out_node.name
                    self.e_reduce[ (v, w) ] = True


            # need to remove the node from the graph rather than just mark the edges are "reduced"?


    def mark_spur_edge(self):

        removed_edges = set()
        for  v in self.nodes:
            if len(self.nodes[v].out_edges) > 1:
                for out_edge in self.nodes[v].out_edges:
                    w = out_edge.out_node.name
                    
                    if len(self.nodes[w].out_edges) == 0 and self.e_reduce[ (v, w) ] != True:
                        self.e_reduce[(v, w)] = True
                        removed_edges.add( (v, w) )
                        v2, w2 = reverse_end(w), reverse_end(v)
                        self.e_reduce[(v2, w2)] = True
                        removed_edges.add( (v2, w2) )

            if len(self.nodes[v].in_edges) > 1:
                for in_edge in self.nodes[v].in_edges:
                    w = in_edge.in_node.name
                    if len(self.nodes[w].in_edges) == 0 and self.e_reduce[ (w, v) ] != True:
                        self.e_reduce[(w, v)] = True
                        removed_edges.add( (w, v) )
                        v2, w2 = reverse_end(w), reverse_end(v)
                        self.e_reduce[(w2, v2)] = True
                        removed_edges.add( (w2, v2) )
        return removed_edges

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
                    v_name, w_name = reverse_end(w.name), reverse_end(v.name)
                    e_reduce[(v_name, w_name)] = True
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
                    v2, w2 = reverse_end(w), reverse_end(v)
                    self.e_reduce[(v2, w2)] = True
                
    def resolve_repeat_edges(self):

        #nxsg = nx.DiGraph()
        #for v, w in self.edges:
        #    if self.e_reduce[(v, w)] != True:
        #        nxsg.add_edge(v, w)
        #nxsg_r = nxsg.reverse()

        edges_to_reduce = []
        nodes_to_test = set()
        for v_n, v in self.nodes.items():
            
            out_nodes = []
            for e in v.out_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    out_nodes.append( e.out_node.name )

            in_nodes = []
            for e in v.in_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    in_nodes.append( e.in_node.name )

            if len(out_nodes) == 1 and len(in_nodes)  == 1:
                nodes_to_test.add(v_n)
        
        for v_n in list( nodes_to_test ):
            
            v = self.nodes[v_n]

            out_nodes = []
            for e in v.out_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    out_nodes.append( e.out_node.name )

            in_nodes = []
            for e in v.in_edges:
                if self.e_reduce[(e.in_node.name, e.out_node.name)] == False:
                    in_nodes.append( e.in_node.name )

            in_node_name = in_nodes[0] 

            for out_edge in self.nodes[in_node_name].out_edges:
                vv = out_edge.in_node.name
                ww = out_edge.out_node.name

                ww_out = self.nodes[ww].out_edges
                v_out = self.nodes[v_n].out_edges
                ww_out_nodes = set( [ n.out_node.name for n in ww_out] )
                v_out_nodes = set(  [ n.out_node.name for n in v_out] )
                o_overlap = len( ww_out_nodes & v_out_nodes )

                ww_in_count = 0
                for e in self.nodes[ww].in_edges:
                    if self.e_reduce[ ( e.in_node.name, e.out_node.name ) ] == False:
                        ww_in_count += 1
                #print ww, ww_in_count 
                #G1 = nx.ego_graph( nxsg,  ww, 3, undirected=False)
                #G2 = nx.ego_graph( nxsg,  v_n, 3, undirected=False)
                #o_overlap = len( set(G1.nodes()) & set(G2.nodes()) )

                if ww != v_n and\
                   self.e_reduce[ (vv, ww) ] == False and\
                   ww_in_count > 1 and\
                   ww not in nodes_to_test and\
                   o_overlap == 0:
                    edges_to_reduce.append( (vv, ww) )

            out_node_name = out_nodes[0]

            for in_edge in self.nodes[out_node_name].in_edges:
                vv = in_edge.in_node.name
                ww = in_edge.out_node.name

                vv_in = self.nodes[vv].in_edges
                v_in = self.nodes[v_n].in_edges
                vv_in_nodes = set( [ n.in_node.name for n in vv_in] )
                v_in_nodes = set(  [ n.in_node.name for n in v_in] )
                i_overlap = len( vv_in_nodes & v_in_nodes )

                vv_out_count = 0
                for e in self.nodes[vv].out_edges:
                    if self.e_reduce[ ( e.in_node.name, e.out_node.name )] == False:
                        vv_out_count += 1
                #print vv, vv_out_count 
                #G1 = nx.ego_graph( nxsg_r,  vv, 3, undirected=False)
                #G2 = nx.ego_graph( nxsg_r,  v_n, 3, undirected=False)
                #i_overlap = len( set(G1.nodes()) & set(G2.nodes()) )

                if vv != v_n and\
                   self.e_reduce[ (vv, ww) ] == False and\
                   vv_out_count > 1 and\
                   vv not in nodes_to_test and\
                   i_overlap == 0:
                    edges_to_reduce.append( (vv, ww) )

        removed_edges = set()
        for e in edges_to_reduce:
            self.e_reduce[e] = True
            removed_edges.add(e)

        return removed_edges

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

def reverse_end( node_id ):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end

def reverse_edge( e ):
    e1, e2 = e
    return reverse_end(e2), reverse_end(e1)

def reverse_path( p ):
    p = p[::-1]
    return [reverse_end(n) for n in p]

    
def find_bundle(sg, edge_data, start_node, depth_cutoff, width_cutoff, length_cutoff):

    tips = set()
    bundle_edges = set()
    bundle_nodes = set()

    local_graph = nx.ego_graph(sg, start_node, depth_cutoff, undirected=False)
    length_to_node = {start_node:0}
    score_to_node = {start_node:0}

    v = start_node
    end_node = start_node

    bundle_nodes.add(v)
    for v, w in local_graph.out_edges(v):
        max_score = 0
        max_length = 0
        if w not in bundle_nodes and reverse_end(w) not in bundle_nodes:
            bundle_edges.add( (v, w) )
            tips.add(w)

    for v in list(tips):
        max_score_node = None
        max_score = 0
        for u, v in local_graph.in_edges(v):
            if u not in length_to_node:
                continue
            score = edge_data[ (u, v) ][4]
            if score > max_score:
                max_score = score
                max_score_node = u

        length_to_node[v] = length_to_node[max_score_node] +  edge_data[ (max_score_node, v) ][3]
        score_to_node[v] = score_to_node[max_score_node] +  edge_data[ (max_score_node, v) ][4]

    for v in list(tips):
        bundle_nodes.add(v)

    depth = 0
    width = 1.0
    converage = False

    #print "start", start_node

    while 1:
        
        #print tips

        depth += 1
        width = 1.0 * len(bundle_edges) / depth

        if depth > 10 and width > width_cutoff:
            converage = False
            break

        if depth > depth_cutoff:
            converage = False
            break
        
        tips_list = list(tips)

        tip_updated = False
        loop_detect = False

        for v in tips_list:

            if len(local_graph.out_edges(v)) == 0: # dead end route
                tip_updated = False
                break

            max_score_node = None
            max_score = 0

            extend_tip = True
            for u, v in local_graph.in_edges(v):


                if u not in length_to_node:
                    extend_tip = False
                    break

                score = edge_data[ (u, v) ][4]

                if score > max_score:

                    max_score = score
                    max_score_node = u
            
            if extend_tip:
            
                length_to_node[v] = length_to_node[max_score_node] +  edge_data[ (max_score_node, v) ][3]
                score_to_node[v] = score_to_node[max_score_node] +  edge_data[ (max_score_node, v) ][4]

                if length_to_node[v] > length_cutoff:
                    converage = False
                    break

                for v, w in local_graph.out_edges(v):
                    if w in length_to_node:
                        loop_detect = True

                    if w not in bundle_nodes and reverse_end(w) not in bundle_nodes:
                        tips.add(w)
                        bundle_edges.add( (v, w) )
                        tip_updated = True

                    if w in bundle_nodes and (v, w) not in bundle_edges:
                        bundle_edges.add( (v, w) )
                        tip_updated = True
                
                tips.remove(v)

            if len(tips) == 1:
                break

        if loop_detect:
            converage = False
            break

        if not tip_updated:
            converage = False
            break

        for v in list(tips):
            bundle_nodes.add(v)

        if len(tips) == 1:

            end_node = tips.pop()
            v = end_node
            max_score_node = None
            max_score = 0

            for u, v in local_graph.in_edges(v):

                if u not in length_to_node:
                    continue

                score = edge_data[ (u, v) ][4]
                if score > max_score:
                    max_score = score
                    max_score_node = u

            if max_score_node == None:
                converage = False
                break

            length_to_node[v] = length_to_node[max_score_node] +  edge_data[ (max_score_node, v) ][3]
            score_to_node[v] = score_to_node[max_score_node] +  edge_data[ (max_score_node, v) ][4]
            converage = True

            break

    data = start_node, end_node, bundle_edges, length_to_node[end_node], score_to_node[end_node], depth
        
    start_node_r, end_node_r = reverse_end(end_node), reverse_end(start_node)
    
    bundle_edge_r = set()
    for v, w in list(bundle_edges):
        vv, ww = reverse_end(w), reverse_end(v)
        bundle_edge_r.add( (vv, ww) )

    data_r = start_node_r, end_node_r, bundle_edge_r, length_to_node[end_node], score_to_node[end_node], depth

    #print converage, data, data_r
    return converage, data, data_r

def generate_string_graph(args):

    overlap_file = args.overlap_file

    contained_reads = set()
    chimer_ids = set()

    filter_reads = False
    
    seqs = set()

    if args.pread_names != "":
        filter_reads = True
        with open(args.pread_names) as f:
            for l in f:
                l = l.strip()
                seqs.add(l)


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
            #if len(l) != 13:
            #    continue
            
            f_id, g_id, score, identity = l[:4]

            if f_id == g_id:  # don't need self-self overlapping
                continue
            
            if filter_reads:

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
            """
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
            """

            overlap_data.append( (f_id, g_id, score, identity,
                                  f_strain, f_start, f_end, f_len,
                                  g_strain, g_start, g_end, g_len) )

            overlap_count[f_id] = overlap_count.get(f_id,0)+1
            overlap_count[g_id] = overlap_count.get(g_id,0)+1
            
    #print "###", len(overlap_data), len(contained_reads)
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
        if f_b > 1:
            if g_b < g_e:
                """
                     f.B         f.E
                  f  ----------->
                  g         ------------->
                            g.B           g.E
                """
                if f_b == 0 or g_e - g_l == 0:
                    continue
                sg.add_edge( "%s:B" % g_id, "%s:B" % f_id, label = (f_id, f_b, 0), 
                                                           length = abs(f_b-0),
                                                           score = -score, 
                                                           identity = identity )
                sg.add_edge( "%s:E" % f_id, "%s:E" % g_id, label = (g_id, g_e, g_l), 
                                                           length = abs(g_e-g_l),
                                                           score = -score,
                                                           identity = identity)
            else:
                """
                     f.B         f.E
                  f  ----------->
                  g         <-------------
                            g.E           g.B           
                """
                if f_b == 0 or g_e == 0:
                    continue
                sg.add_edge( "%s:E" % g_id, "%s:B" % f_id, label = (f_id, f_b, 0), 
                                                           length = abs(f_b -0),
                                                           score = -score,
                                                           identity = identity)
                sg.add_edge( "%s:E" % f_id, "%s:B" % g_id, label = (g_id, g_e, 0), 
                                                           length = abs(g_e- 0),
                                                           score = -score,
                                                           identity = identity)
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
                sg.add_edge( "%s:B" % f_id, "%s:B" % g_id, label = (g_id, g_b, 0), 
                                                           length = abs(g_b - 0),
                                                           score = -score,
                                                           identity = identity)
                sg.add_edge( "%s:E" % g_id, "%s:E" % f_id, label = (f_id, f_e, f_l), 
                                                           length = abs(f_e-f_l),
                                                           score = -score,
                                                           identity = identity)
            else:
                """
                                    f.B         f.E
                  f                 ----------->
                  g         <-------------
                            g.E           g.B           
                """
                if g_b - g_l == 0 or f_e - f_l ==0:
                    continue
                sg.add_edge( "%s:B" % f_id, "%s:E" % g_id, label = (g_id, g_b, g_l), 
                                                           length = abs(g_b - g_l),
                                                           score = -score,
                                                           identity = identity)
                sg.add_edge( "%s:B" % g_id, "%s:E" % f_id, label = (f_id, f_e, f_l), 
                                                           length = abs(f_e - f_l),
                                                           score = -score,
                                                           identity = identity)


    sg.init_reduce_dict()
    #if not args.disable_chimer_prediction:
    #    sg.mark_chimer_edges()
    #sg.mark_spur_edge()
    

    sg.mark_tr_edges() # mark those edges that transitive redundant

    if DEBUG_LOG_LEVEL > 1:
        print sum( [1 for c in sg.e_reduce.values() if c == True] )
        print sum( [1 for c in sg.e_reduce.values() if c == False] )


    #sg.mark_best_overlap() # mark those edges that are best overlap edges
    removed_edges = sg.resolve_repeat_edges()  
    #removed_edges = set()

    spur_edges = sg.mark_spur_edge()

    if DEBUG_LOG_LEVEL > 1:
        print sum( [1 for c in sg.e_reduce.values() if c == False] )

    #max_score = max([ sg.edges[ e ].attr["score"] for e in sg.edges if sg.e_reduce[e] != True ])
    out_f = open("sg_edges_list", "w")
    nxsg = nx.DiGraph()
    edge_data = {}
    for v, w in sg.edges:
        if sg.e_reduce[(v, w)] != True or (v, w) in removed_edges or (v, w) in spur_edges:
            e = sg.edges[ (v, w) ]
            rid, sp, tp = e.attr["label"]
            score = e.attr["score"]
            identity = e.attr["identity"]
            length = abs(sp-tp)


            if  sg.e_reduce[(v, w)] != True:
                type_ = "G"
                label = "%s:%d-%d" % (rid, sp, tp)
                nxsg.add_edge(v, w, label = label, length = length, score = score)

            elif (v, w) in removed_edges:
                type_ = "R"

            elif (v, w) in spur_edges:
                type_ = "S"

            print >>out_f, v, w, rid, sp, tp, score, identity, type_
            edge_data[ (v, w) ] = (rid, sp, tp, length, score, identity, type_)

        
    out_f.close()
    nxsg_r = nxsg.reverse()    

    return nxsg, nxsg_r, edge_data



def construct_compound_paths(sg, sg_r, edge_data):

    source_nodes = set()
    sink_nodes = set()
    simple_nodes = set()
    branch_nodes = set()

    all_nodes = sg.nodes()
    for n in all_nodes:
        in_degree = len( sg.in_edges(n) )
        out_degree = len( sg.out_edges(n) )
        if in_degree == 0:
            source_nodes.add(n)
        if out_degree == 0:
            sink_nodes.add(n)
        if in_degree == 1 and out_degree == 1:
            simple_nodes.add(n)
        if in_degree > 1 or out_degree > 1:
            branch_nodes.add(n)

    #print "#", len(all_nodes),len(source_nodes), len(sink_nodes), len(simple_nodes), len(branch_nodes)
    compound_paths_ = {}

    compound_path_start_at = {}
    compound_path_end_at = {}
    for p in list(branch_nodes):
        if sg.out_degree(p) > 1:
            coverage, data, data_r =  find_bundle(sg, edge_data, p, 128, 8, 500000)
            if coverage == True:
                start_node, end_node, bundle_edges, length, score, depth = data
                compound_paths_[ (start_node, "NA", end_node) ] = ( 1.0*len(bundle_edges)/depth, length, score, bundle_edges )
               
                compound_path_start_at.setdefault( start_node, [] )
                compound_path_start_at[start_node].append( ( len(bundle_edges), (start_node, "NA", end_node) ) )
                
                compound_path_end_at.setdefault( end_node, [] )
                compound_path_end_at[end_node].append( ( len(bundle_edges), (start_node, "NA", end_node) ) )
                
                start_node, end_node, bundle_edges, length, score, depth = data_r
                compound_paths_[ (start_node, "NA", end_node) ] = ( 1.0*len(bundle_edges)/depth, length, score, bundle_edges )

                compound_path_start_at.setdefault( start_node, [] )
                compound_path_start_at[start_node].append( ( len(bundle_edges), (start_node, "NA", end_node) ) )

                compound_path_end_at.setdefault( end_node, [] )
                compound_path_end_at[end_node].append( ( len(bundle_edges), (start_node, "NA", end_node) ) )

    duplicate_edges = set()
    for start_node in compound_path_start_at:
        if len(compound_path_start_at[start_node]) == 1:
            continue

        compound_path_start_at[start_node].sort()
        
        max_path_id = compound_path_start_at[start_node][-1][1]

        max_cover_edges = set(compound_paths_[max_path_id][3])

        for l_e, path_id in compound_path_start_at[start_node][:-1]:
            if path_id == max_path_id:
                continue
            edges = set(compound_paths_[path_id][3])
            if len( max_cover_edges | edges ) == len( max_cover_edges ) : #edges is a subset if max_cover_edges
                #print "remove", path_id, "covered by", compound_path_start_at[start_node][-1][1]
                duplicate_edges.add( path_id )

    for end_node in compound_path_end_at:
        if len(compound_path_end_at[end_node]) == 1:
            continue

        compound_path_end_at[end_node].sort()
        
        max_path_id = compound_path_end_at[end_node][-1][1]

        max_cover_edges = set(compound_paths_[max_path_id][3])

        for l_e, path_id in compound_path_end_at[end_node][:-1]:
            if path_id == max_path_id:
                continue
            edges = set(compound_paths_[path_id][3])
            if len( max_cover_edges | edges ) == len( max_cover_edges ) : #edges is a subset if max_cover_edges
                #print "remove", path_id, "covered by", compound_path_end_at[end_node][-1][1]
                duplicate_edges.add( path_id )

    #for path_id in list(duplicate_edges):
    #    print "remove",path_id
    compound_paths = {}
    for k, v in compound_paths_.items():
        if k in duplicate_edges:
            continue
        compound_paths[ k ] = v

    return compound_paths

if __name__ == "__main__":

    import argparse
    
    parser = argparse.ArgumentParser(description='a example string graph assembler that is desinged for handling diploid genomes')
    parser.add_argument('overlap_file', help='a file that contains the overlap information.')
    parser.add_argument('--pread_names', default = "", help='the file that contains the name of the sequence to be assembled')

    parser.add_argument('--min_len', type=int, default=4000, 
                        help='minimum length of the reads to be considered for assembling')
    parser.add_argument('--min_idt', type=float, default=96,
                        help='minimum alignment identity of the reads to be considered for assembling')
    parser.add_argument('--disable_chimer_prediction', action="store_true", default=False,
                        help='you may want to disable this as some reads can be falsely identified as chimers in low coverage case')

    args = parser.parse_args()


    # transitivity reduction, remove spurs, remove putative edges caused by repeats
    sg, sg_r, edge_data = generate_string_graph(args)

    # identify compound_paths
    compound_paths = construct_compound_paths(sg, sg_r, edge_data)

    simple_paths = {}
    dual_path = {}

    for s, v, t in compound_paths:
        dual_path[ (s, v, t)  ] = (reverse_end(t), v, reverse_end(s))
        dual_path[ (reverse_end(t), v, reverse_end(s)) ] = (s, v, t)
        assert (reverse_end(t), v, reverse_end(s)) in compound_paths


    masked_edges = set()
    for e in compound_paths:
        edges = compound_paths[e][3]
        for e in edges:
            masked_edges.add( (e[0], e[1]) )
            #masked_edges.add( (reverse_end(e[1]), reverse_end( e[0]) ) )
        
    for s,t in list(masked_edges):
        assert (reverse_end(t), reverse_end(s)) in masked_edges


    # create a string graph without compound edges
    sg2 = nx.DiGraph()

    for v, w in edge_data:

        assert (reverse_end(w), reverse_end(v)) in edge_data
        
        if (v, w) in masked_edges:
            continue

        rid, sp, tp, length, score, identity, type_ = edge_data[ (v, w) ]
        if type_ != "G":
            continue

        label = "%s:%d-%d" % (rid, sp, tp)
        sg2.add_edge( v, w, label = label, length = length, score = score)

    for s, v, t in compound_paths:
        width, length, score, bundle_edges =  compound_paths[ (s, v, t) ] 
        sg2.add_edge( s, t, label = "c:%s~%s" % (s, t), length = length, score = score)
                

    # utg construction phase 1, identify all simple paths
    s_nodes = set()
    t_nodes = set()
    simple_nodes = set()

    all_nodes = sg2.nodes()
    for n in all_nodes:
        in_degree = len( sg2.in_edges(n) )
        out_degree = len( sg2.out_edges(n) )
        if in_degree == 1 and out_degree == 1:
            simple_nodes.add(n)
        else:
            if out_degree != 0:
                s_nodes.add(n)
            if in_degree != 0:
                t_nodes.add(n)

    free_edges = set(sg2.edges())
    
    for s, v, t in compound_paths:
        free_edges.remove( (s, t) )

    for v,w in free_edges:
        if (reverse_end(w), reverse_end(v) ) not in free_edges:
            print v,w
            print oreverse_end(w), reverse_end(v)

    while len(free_edges) != 0:
        if len(s_nodes) != 0:
            n = s_nodes.pop()
        else:
            e = free_edges.pop()
            free_edges.add(e)
            n = e[0]


        path = []
        path_length =0
        path_score = 0 
        for v, w in sg2.out_edges(n):
            if (v, w) not in free_edges:
                continue
            rv = reverse_end(v)
            rw = reverse_end(w)
            #if (rw, rv) not in free_edges:
            #    continue

            path_length = 0
            path_score = 0
            v0 = v
            w0 = w
            path = [v, w]
            path_edges = set()
            path_edges.add( (v, w) )
            path_length += edge_data[ (v, w) ][3]
            path_score += edge_data[ (v, w) ][4]
            free_edges.remove( (v, w) )

            r_path_length = 0
            r_path_score = 0
            rv0 = rv
            rw0 = rw
            r_path = [rv, rw] # need to reverse again
            r_path_edges = set()
            r_path_edges.add( (rw, rv) )
            r_path_length += edge_data[ (rw, rv) ][3]
            r_path_score += edge_data[ (rw, rv) ][4]
            free_edges.remove( (rw, rv) )

            while w in simple_nodes:
                w, w_ = sg2.out_edges(w)[0]
                if (w, w_) not in free_edges:
                    break
                rw_, rw = reverse_end(w_), reverse_end(w)

                if ( rw_, rw ) in path_edges:
                    break

                path.append(w_)
                path_edges.add( (w, w_) )
                path_length += edge_data[ (w, w_) ][3]
                path_score += edge_data[ (w, w_) ][4]
                free_edges.remove( (w, w_) )
                
                r_path.append(rw_)
                r_path_edges.add( (rw_, rw) )
                r_path_length += edge_data[ (rw_, rw) ][3]
                r_path_score += edge_data[ (rw_, rw) ][4]
                free_edges.remove( (rw_, rw) )
                

                w = w_

            simple_paths[ (v0, w0, path[-1]) ] = path_length, path_score, path
            r_path.reverse()
            assert r_path[0] == reverse_end(path[-1])
            simple_paths[ (r_path[0], rw0, rv0) ] = r_path_length, r_path_score, r_path

            dual_path[ (r_path[0], rw0, rv0) ] = (v0, w0, path[-1])
            dual_path[ (v0, w0, path[-1]) ] = (r_path[0], rw0, rv0)


    #phase 2

    #convert simple paths that have the same source and sink to a compound path
    
    st_path = {}
    for s, v, t in simple_paths:
        length, score, path = simple_paths[ (s, v, t) ]
        st_path.setdefault( (s,t), [] )
        st_path[ (s, t) ].append( ( (s, v, t), length, score,  path ) )
    
    for s, v, t in compound_paths:
        width, length, score, edges = compound_paths[ (s, v, t) ]
        st_path.setdefault( (s,t), [] )
        st_path[ (s, t) ].append( ( (s, v, t), length, score,  edges ) )

    simple_paths = {}
    compound_paths = {}

    for s, t in st_path:
        if len( st_path[ (s, t) ] ) == 1:
            svt, length, score, path_or_edges = st_path[ (s, t) ][0]
            if svt[1] == "NA":
                compound_paths[ svt ] = length, score, path_or_edges
            else:
                simple_paths[ svt ] = length, score, path_or_edges
        else:
            edges = []
            total_score = 0
            max_length = 0
            for svt, length, score, path_or_edges in st_path[ (s, t) ]:
                ss, vv, tt = svt
                total_score += score
                if length > max_length:
                    max_length = length
                if vv[1] != "NA": #a path
                    path = path_or_edges
                    n1 = path[0]
                    for i in xrange(1, len(path_or_edges)):
                        n2 = path[i]
                        edges.append( (n1, n2) )
                        n1 = n2
                else:
                    edges.extend( path_or_edges )
            compound_paths[ (ss, "NA", tt) ] = max_length, total_score, edges 
            dual_path[ (v0, w0, path[-1]) ] = (r_path[0], rw0, rv0)

    for s, v, t in compound_paths:
        dual_path[ (s, v, t)  ] = (reverse_end(t), v, reverse_end(s))
        dual_path[ (reverse_end(t), v, reverse_end(s)) ] = (s, v, t)

    ug = nx.MultiDiGraph()
    u_edge_data = {}
    circular_path = set()

    for s, v, t in simple_paths:
        length, score, path = simple_paths[ (s, v, t) ]
        u_edge_data[ (s, t, v) ] = (length, score, path)
        if s != t:
            ug.add_edge(s, t, key = v, type_ = "simple", via = v, length = length, score = score)
        else:
            circular_path.add( (s, t, v) )

    for s, v, t in compound_paths:
        length, score, edges = compound_paths[ (s, v, t) ]
        u_edge_data[ (s, t, v) ] = (length, score, edges)
        assert v == "NA"
        if s != t:
            ug.add_edge(s, t, key = v, type_ = "compound", via = v, length = length, score = score)
        else:
            circular_path.add( (s, t, v) )


    # identify spurs in the utg graph
    # Currently, we use ad-hoc logic filtering out shorter utg, but we ca
    # add proper alignment comparison later to remove redundant utgs 
    utg_spurs = set()

    all_nodes = ug.nodes()
    for n in all_nodes:
        in_degree = len( ug.in_edges(n) )
        out_degree = len( ug.out_edges(n) )
        if in_degree == 1 and out_degree == 1:
            simple_nodes.add(n)
        else:
            if out_degree != 0:
                s_nodes.add(n)
            if in_degree != 0:
                t_nodes.add(n)
    
    for n in s_nodes:
        max_length = 0
        more_connection = 0
        for s, t, v in ug.out_edges(n, keys=True):
            length, score, edges = u_edge_data[ (s, t, v) ]
            if length > max_length:
                max_length = length
            if ug.out_degree(t) > 0:
                more_connection = 1

        for s, t, v in ug.out_edges(n, keys=True):
            length, score, edges = u_edge_data[ (s, t, v) ]
            if length < max_length * 0.25 and length < 50000:  
                utg_spurs.add( ( s, t, v ) )
                ds, dv, dt = dual_path[ (s, v, t) ]
                utg_spurs.add( (ds, dt, dv) ) 
            elif length < 50000 and more_connection == 1 and ug.out_degree(t) == 0:
                utg_spurs.add( ( s, t, v ) )
                ds, dv, dt = dual_path[ (s, v, t) ]
                utg_spurs.add( (ds, dt, dv) ) 


    for n in t_nodes:
        max_length = 0
        more_connection = 0
        for s, t, v in ug.in_edges(n, keys=True):
            length, score, edges = u_edge_data[ (s, t, v) ]
            if length > max_length:
                max_length = length
            if ug.in_degree(s) > 0:
                more_connection = 1
        for s, t, v in ug.in_edges(n, keys=True):
            length, score, edges = u_edge_data[ (s, t, v) ]
            if length < max_length * 0.25 and length < 50000:
                utg_spurs.add( ( s, t, v ) )
                ds, dv, dt = dual_path[ (s, v, t) ]
                utg_spurs.add( (ds, dt, dv) ) 
            elif length < 50000 and more_connection == 1 and ug.in_degree(s) == 0:
                utg_spurs.add( ( s, t, v ) )
                ds, dv, dt = dual_path[ (s, v, t) ]
                utg_spurs.add( (ds, dt, dv)  ) 



    with open("utg_data","w") as f:
        for s, t, v in u_edge_data:
            if (s, t, v) in utg_spurs:
                is_spur = "+"
                try:
                    ug.remove_edge(s, t, key=v)
                except: #some self-dual edge can be already removed
                    pass
            else:
                is_spur = "-"

            length, score, path_or_edges = u_edge_data[ (s, t, v) ]
            if v == "NA":
                type_ = "edges"
                path_or_edges = ".".join( [ vv+"~"+ww for vv, ww in path_or_edges ] )
            else:
                type_ = "path"
                path_or_edges = "~".join( path_or_edges )
            print >>f, s, v, t, type_, is_spur, length, score, path_or_edges



    # contig construction from utgs

    s_nodes = set()
    t_nodes = set()
    simple_nodes = set()

    all_nodes = ug.nodes()
    for n in all_nodes:
        in_degree = len( ug.in_edges(n) )
        out_degree = len( ug.out_edges(n) )
        if in_degree == 1 and out_degree == 1:
            simple_nodes.add(n)
        else:
            if out_degree != 0:
                s_nodes.add(n)
            if in_degree != 0:
                t_nodes.add(n)

    free_edges = set()
    for s, t, v, d in ug.edges(keys=True, data=True):
        free_edges.add( (s, t, v) )

    ctg_id = 0

    ctg_paths = open("ctg_paths","w")

    while len(free_edges) != 0:
        if len(s_nodes) != 0:
            n = s_nodes.pop()
        else:
            e = free_edges.pop()
            free_edges.add(e)
            n = e[0]

        path = []
        path_length =0
        path_score = 0
        r_path = []
        r_path_length =0
        r_path_score = 0
        for s, t, v, d in ug.out_edges(n, keys=True, data=True):
            #v = d["via"]
            if (s, t, v) not in free_edges:
                continue

            rs = reverse_end(s)
            rt = reverse_end(t)
            if v != "NA":
                rv = reverse_end(v)
            else:
                rv = v
            
            path_length = 0
            path_score = 0
            s0, t0, v0 = s, t, v
            the_edge = (s, t, v) 
            path = [ the_edge ]
            path_edges = set()
            path_edges.add( the_edge )
            path_length += u_edge_data[ the_edge ][0]
            path_score += u_edge_data[ the_edge ][1]
            free_edges.remove( the_edge )
            
            r_path_length = 0
            r_path_score = 0
            rs0, rt0, rv0 = rs, rt, rv
            the_edge = (rt, rs, rv)
            r_path = [ the_edge ]
            r_path_edges = set()
            r_path_edges.add( the_edge )
            r_path_length += u_edge_data[ the_edge ][0]
            r_path_score += u_edge_data[ the_edge ][1]
            try:
                free_edges.remove( (rt, rs, rv) )
            except:
                pass

            while t in simple_nodes:

                t, t_, d = ug.out_edges( t, data = True )[0]
                v = d["via"]

                if (t, t_, v) not in free_edges:
                    break

                if ( reverse_end(t_), reverse_end(t) ) in path_edges:
                    break

                the_edge = (t, t_, v)
                path.append( the_edge )
                path_edges.add( (t, t_ ) )
                path_length += u_edge_data[ the_edge ][0]
                path_score += u_edge_data[ the_edge ][1]
                #print path_length
                free_edges.remove( the_edge )
                if v != "NA": 
                    rt, rt_, rv = reverse_end(t), reverse_end(t_), reverse_end(v)
                else:
                    rt, rt_, rv = reverse_end(t), reverse_end(t_), "NA"

                the_edge = (rt_, rt, rv)
                r_path.append( the_edge )
                r_path_edges.add( (rt_, rt) )
                r_path_length += u_edge_data[ the_edge ][0]
                r_path_score += u_edge_data[ the_edge ][1]
                try:
                    free_edges.remove( the_edge ) 
                except:
                    pass

                t = t_

            print >> ctg_paths, "%06dF" % ctg_id, "ctg_linear", s0+"~"+v0+"~"+t0, path[-1][1], path_length, path_score, ".".join([ c[0]+"~"+c[2]+"~"+c[1] for c in path ] )
            r_path.reverse()
            print >> ctg_paths, "%06dR" % ctg_id, "ctg_linear", r_path[0][0]+"~"+r_path[0][2]+"~"+r_path[0][1], r_path[-1][1], r_path_length, r_path_score, ".".join([ c[0]+"~"+c[2]+"~"+c[1] for c in r_path ] )
            ctg_id += 1

        


    for s, t, v in list(circular_path):
        length, score, path = u_edge_data[ (s, t, v) ]
        print >> ctg_paths, "%6d" % ctg_id, "ctg_circular", s+"~"+v+"~"+t, t, length, score, s+"~"+v+"~"+t
        ctg_id += 1

    ctg_paths.close()

