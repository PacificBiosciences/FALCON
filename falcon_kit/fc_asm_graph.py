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
from FastaReader import FastaReader

RCMAP = dict(zip("ACGTacgtNn-","TGCAtgcaNn-"))

def reverse_end( node_id ):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end

class AsmGraph(object):

    def __init__(self, sg_file, utg_file, ctg_file):
        self.sg_edges = {}
        self.sg_edge_seqs = {}
        self.utg_data = {}
        self.ctg_data ={}
        self.utg_to_ctg = {}
        self.node_to_ctg = {}
        self.node_to_utg = {}

        self.load_sg_data(sg_file)
        self.load_utg_data(utg_file)
        self.load_ctg_data(ctg_file)

        self.build_node_map()

    def load_sg_data(self, sg_file):

        with open(sg_file) as f:
            for l in f:
                l = l.strip().split()
                v, w = l[0:2]
                seq_id, b, e = l[2:5]
                b, e = int(b), int(e)
                score, idt = l[5:7]
                score, idt = int(score), float(idt)
                type_ = l[7]
                self.sg_edges[ (v, w) ] = ( (seq_id, b, e), score, idt, type_)

    def load_sg_seq(self, fasta_fn):

        all_read_ids = set() # read ids in the graph

        for v, w in self.sg_edges:
            type_ = self.sg_edges[ (v, w) ][-1]
            if type_ != "G":
                continue
            v = v.split(":")[0]
            w = w.split(":")[0]
            all_read_ids.add(v)
            all_read_ids.add(w)

        seqs = {}
        # load all p-read name into memory
        f = FastaReader(fasta_fn)
        for r in f:
            if r.name not in all_read_ids:
                continue
            seqs[r.name] = r.sequence.upper()


        for v, w in self.sg_edges:
            seq_id, s, t = self.sg_edges[ (v, w) ][0]
            type_ = self.sg_edges[ (v, w) ][-1]

            if type_ != "G":
                continue

            if s < t:
                e_seq = seqs[ seq_id ][ s:t ]
            else:
                e_seq = "".join([ RCMAP[c] for c in seqs[ seq_id ][ s:t:-1 ] ])
            self.sg_edge_seqs[ (v, w) ] = e_seq

    def get_seq_from_path(self, path):
        if len(self.sg_edge_seqs) == 0:
            return ""
        v = path[0]
        seqs = []
        for w in path[1:]:
            seqs.append( self.sg_edge_seqs[ (v, w) ] )
            v = w
        return "".join(seqs)


    def load_utg_data(self, utg_file):

        with open(utg_file) as f:
            for l in f:
                l = l.strip().split()
                s, v, t = l[0:3]
                type_, length, score = l[3:6]
                length, score = int(length), int(score)
                path_or_edges = l[6]
                self.utg_data[ (s,t,v) ] = ( type_, length, score, path_or_edges)


    def load_ctg_data(self, ctg_file):

        with open(ctg_file) as f:
            for l in f:
                l = l.strip().split()
                ctg_id, ctg_type = l[0:2]
                start_edge = l[2]
                end_node = l[3]
                length = int(l[4])
                score = int(l[5])
                path = tuple( ( e.split("~") for e in l[6].split("|") ) )
                self.ctg_data[ ctg_id ] = ( ctg_type, start_edge, end_node,  length, score, path )
                for u in path:
                    s, v, t = u
                    #rint s,v,t
                    type_, length, score, path_or_edges =  self.utg_data[ (s,t,v) ]
                    if type_ != "compound":
                        self.utg_to_ctg[ (s, t, v) ] = ctg_id
                    else:
                        for svt in path_or_edges.split("|"):
                            s, v, t = svt.split("~")
                            self.utg_to_ctg[ (s, t, v) ] = ctg_id


    def get_sg_for_utg(self, utg_id):
        sg = nx.DiGraph()
        type_, length, score, path_or_edges =  self.utg_data[ utg_id ]
        if type_ == "compound":
            for svt in path_or_edges.split("|"):
                s, v, t = svt.split("~")
                type_, length, score, one_path =  self.utg_data[ (s, t, v) ]
                one_path = one_path.split("~")
                sg.add_path(one_path)
        else:
            one_path = path_or_edges.split("~")
            sg.add_path(one_path)
        return sg


    def get_sg_for_ctg(self, ctg_id):
        sg = nx.DiGraph()
        utgs = []
        path = self.ctg_data[ctg_id][-1]
        for s, v, t in path:
            type_, length, score, path_or_edges =  self.utg_data[ (s, t, v) ]
            utgs.append( (type_, path_or_edges) )

        for t, utg in utgs:
            if t == "simple":
                one_path = utg.split("~")
                sg.add_path(one_path)
            elif t == "compound":
                for svt in utg.split("|"):
                    s, v, t = svt.split("~")
                    type_, length, score, one_path =  self.utg_data[ (s, t, v) ]
                    one_path = one_path.split("~")
                    sg.add_path(one_path)

        return sg


    def build_node_map(self):

        for ctg_id in self.ctg_data:
            sg = self.get_sg_for_ctg( ctg_id )
            for n in sg.nodes():
                self.node_to_ctg.setdefault(n, set())
                self.node_to_ctg[n].add(ctg_id)


        for u_id in self.utg_data:
            if self.utg_data[u_id][0] == "compound":
                continue
            sg = self.get_sg_for_utg( u_id )
            for n in sg.nodes():
                self.node_to_utg.setdefault(n, set())
                self.node_to_utg[n].add( u_id )
