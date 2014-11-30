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

from fc_asm_graph import AsmGraph

G_asm = AsmGraph("sg_edges_list", "utg_data", "ctg_paths")

sg_edges = G_asm.sg_edges
node_to_ctg = G_asm.node_to_ctg
node_to_utg = G_asm.node_to_utg

ctg_data = G_asm.ctg_data
utg_data = G_asm.utg_data

ctg_pair_links = {}
for v, w in sg_edges.keys():
    if v in node_to_ctg and w in node_to_ctg:
        for ctg1 in list(node_to_ctg[v]):
            for ctg2 in list(node_to_ctg[w]):
                if ctg1 == ctg2:
                    continue
                ctg_pair_links.setdefault((ctg1, ctg2), set())
                ctg_pair_links[ (ctg1, ctg2) ].add( (v,w) )    

                    
utg_pair_links = {}
for v, w in sg_edges.keys():
    if v in node_to_utg and w in node_to_utg:
        for u1 in list(node_to_utg[v]):
            for u2 in list(node_to_utg[w]):
                if u1 == u2:
                    continue
                utg_pair_links.setdefault((u1, u2), set())
                utg_pair_links[(u1,u2)].add( (v, w) )


for ctg1, ctg2 in ctg_pair_links:
    links = ctg_pair_links[ ( ctg1, ctg2 ) ]
    count = len(links)
    if count > 0:
        path1 = ctg_data[ctg1][-1][-5:]
        path2 = ctg_data[ctg2][-1][:5]
        utg1 = []
        utg2 = []
        for s1, v1, t1 in path1:
            u1 = (s1, t1, v1)
            type_, length, score, path_or_edges =  utg_data[ u1 ]
            if type_ == "compound":
                for u in path_or_edges.split("|"):
                    ss, vv, tt = u.split("~")
                    utg1.append( (ss, tt, vv) )
            else:
               utg1.append(u1)
        for s2, v2, t2 in path2:
            u2 = (s2, t2, v2)
            type_, length, score, path_or_edges =  utg_data[ u2 ]
            if type_ == "compound":
                for u in path_or_edges.split("|"):
                    ss, vv, tt = u.split("~")
                    utg2.append( (ss, tt, vv) )
            else:
               utg2.append(u2) 
        #print path1
        #print path2
        #print len(utg1), len(utg2)
        for u1 in utg1:
            for u2 in utg2:
                u1 = tuple(u1)
                u2 = tuple(u2)
                c = utg_pair_links.get( (u1, u2), set() )
                if len(c) == 0:
                    continue
                s1,t1,v1 = u1
                s2,t2,v2 = u2
                len_1 = ctg_data[ ctg1 ][ 3 ]
                len_2 = ctg_data[ ctg2 ][ 3 ]
                print ctg1, ctg2, len_1, len_2, len(utg1), len(utg2), len(links), "~".join( (s1,v1,t1) ),  "~".join( (s2,v2,t2) ), len(c)
        


