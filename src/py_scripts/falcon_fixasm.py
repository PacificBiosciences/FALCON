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

import networkx as nx
from pbcore.io import FastaReader

def neighbor_bound(G, v, w, radius):
    g1 = nx.ego_graph(G, v, radius=radius, undirected=False)
    g2 = nx.ego_graph(G, w, radius=radius, undirected=False)
    if len(g1) < radius or len(g2) < radius:
        return True
    print v, len(g1), w, len(g2), radius
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
            if neighbor_bound(G, v, w, 20) == False:
                is_branch = True
                break
        if is_branch == True:
            break
    return is_branch


def get_r_path(r_edges, u_path):
    tiling_path = []
    pos = 0
     
    for i in range( len(u_path) - 1): 
        v, w = u_path[i:i+2]
        r_edge_label, overlap = r_edges[ (v, w) ]
        r_edge_seq_id, range_ = r_edge_label.split(":")
        range_ = range_.split("-")
        s, e = int(range_[0]), int(range_[1])
        pos += abs(e-s)
        tiling_path.append( (pos, w, s, e) )
    return tiling_path

def get_seq(u_edges, r_edges, path):
    subseqs = []
    pos = []
    cur_pos = 0
    full_tiling_path = []

    for i in range( len(path) - 1):
        v, w = path[i:i+2]
        pos.append( (v, cur_pos) )
        uedges = u_edges[ (v, w) ]
        uedges.sort( key= lambda x: len(x[0]) )
        subseqs.append( uedges[-1][1] )
        r_path = get_r_path( r_edges, uedges[-1][0].split("-") )
        r_path = [ ( x[0] + cur_pos, x[1], x[2], x[3]) for x in r_path ]
        full_tiling_path.extend( r_path )
        cur_pos += len( uedges[-1][1] )
    pos.append( (w, cur_pos) ) 
    return "".join(subseqs), pos, full_tiling_path


u_edges = {}
with open("unit_edges.dat") as f:
    for l in f:
        v, w, path, seq = l.strip().split()
        u_edges.setdefault( (v, w), [] )
        u_edges[ (v, w) ].append( (path, seq) )
len(u_edges)


r_edges = {}
with open("edges_list") as f:
    for l in f:
        v, w, edge_label, overlap = l.strip().split()
        r_edges[ (v, w) ] = (edge_label, int(overlap) ) 


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
pos_f = open("primary_tigs_node_pos_c", "w")
tiling_path_f = open("all_tiling_path_c", "w")
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
                seq, pos, full_tiling_path = get_seq(u_edges, r_edges, c_path)
                for p, w, s, e in full_tiling_path:
                    print >> tiling_path_f, "%s_%02d" % (name, sub_idx), p, w, s, e
                #if len(full_tiling_path) <= 5:
                #    continue
                print >>out_f, ">%s_%02d" % (name, sub_idx)
                print >>out_f, seq
                print >>path_f, ">%s_%02d" % (name, sub_idx), " ".join(c_path)
                #print c_path
                for node, p in pos:
                    print >> pos_f, "%s_%02d %s %d" % (name, sub_idx, node, p)
                c_path = [v]
                sub_idx += 1
            else:
                c_path.append(v)
                
        if len(c_path) > 1:
            seq, pos, full_tiling_path = get_seq(u_edges, r_edges, c_path)
            for p, w, s, e in full_tiling_path:
                print >> tiling_path_f, "%s_%02d" % (name, sub_idx), p, w, s, e
            if len(full_tiling_path) <= 5:
                continue
            print >>out_f, ">%s_%02d" % (name, sub_idx)
            print >>out_f, seq
            print >>path_f, ">%s_%02d" % (name, sub_idx), " ".join(c_path)
            for node, p in pos:
                print >> pos_f, "%s_%02d %s %d" % (name, sub_idx, node, p)

with open("all_tigs_paths") as f:
    for l in f:
        l = l.strip().split()
        name = l[0][1:]
        name = name.split("-")
        if name[1] == "0000":
            continue
        if len(name) == 2:
            path = l[1:]
            seq, pos, full_tiling_path = get_seq(u_edges, r_edges, path)
            for p, w, s, e in full_tiling_path:
                print >> tiling_path_f, "%s" % ("-".join(name)), p, w, s, e
        else:
            path = l[1:]
            full_tiling_path = get_r_path(r_edges, path)
            for p, w, s, e in full_tiling_path:
                print >> tiling_path_f, "%s" % ("-".join(name)), p, w, s, e

            
path_f.close()
tiling_path_f.close()
pos_f.close()
