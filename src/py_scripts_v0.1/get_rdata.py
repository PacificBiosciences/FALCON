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

import sys
import glob
#import pkg_resources
import uuid
from datetime import datetime

from collections import Counter
from multiprocessing import Pool
#from pbtools.pbdagcon.q_sense import *
import os

"""
try:
    __p4revision__ = "$Revision: #4 $"
    __p4change__ = "$Change: 121571 $"
    revNum = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
    changeNum = int(__p4change__.strip("$").split(":")[-1])
    __version__ = "%s-r%d-c%d" % ( pkg_resources.require("pbtools.pbhgap")[0].version, revNum, changeNum )
except:
    __version__ = "pbtools.hbar-dtk-github"
"""

query_fasta_fn = sys.argv[1]
target_fasta_fn = sys.argv[2]
m4_fofn = sys.argv[3]
bestn = int(sys.argv[4])
group_id = int(sys.argv[5])
num_chunk = int(sys.argv[6])
min_cov = int(sys.argv[7])
max_cov = int(sys.argv[8])
trim_align = int(sys.argv[9])
trim_plr = int(sys.argv[10])


rmap = dict(zip("ACGTNacgt-","TGCANntgca-"))
def rc(seq):
    return "".join([rmap[c] for c in seq[::-1]])

"""0x239fb832/0_590 0x722a1e26 -1843 81.6327 0 62 590 590 0 6417 6974 9822 254 11407 -74.5375 -67.9 1"""
query_to_target = {}
with open(m4_fofn) as fofn:
    for fn in fofn:
        fn = fn.strip()
        with open(fn) as m4_f:
            for l in m4_f:
                d = l.strip().split()
                id1, id2 = d[:2]
                #if -noSplitSubread not used, we will need the following line    
                #id1 = id1.split("/")[0]
                if id1 == id2:
                    continue
                if hash(id2) % num_chunk != group_id:
                    continue
                if int(d[2]) > -1000: continue
                if int(d[11]) < 4000: continue
                query_to_target.setdefault(id1, [])
                query_to_target[id1].append( (int(d[2]), l) )

target_to_query = {}
for id1 in query_to_target:
    query_to_target[id1].sort()
    rank = 0
    for s, ll in query_to_target[id1][:bestn]:
        l = ll.strip()
        d = l.split()
        id1, id2 = d[:2]
        target_to_query.setdefault(id2,[])
        target_to_query[id2].append( ( (int(d[5])-int(d[6]), int(d[2])), l ) )
        #target_to_query[id2].append( ( int(d[2]), l ) )
        #rank += 1

from pbcore.io import FastaIO
query_data = {}
with open(query_fasta_fn) as fofn:
    for fa_fn in fofn:
        fa_fn = fa_fn.strip()
        f_s = FastaIO.FastaReader(fa_fn)
        for s in f_s:
            id1 = s.name
            if id1 not in query_to_target:
                continue
            query_data[id1]=s.sequence
        f_s.file.close()

target_data = {}
with open(target_fasta_fn) as fofn:
    for fa_fn in fofn:
        fa_fn = fa_fn.strip()
        f_s = FastaIO.FastaReader(fa_fn)
        for s in f_s:
            id2 = s.name
            if hash(id2) % num_chunk != group_id:
                continue
            target_data[id2]=s.sequence
        f_s.file.close()


ec_data = []
base_count = Counter()
r_count =0

for id2 in target_to_query:
    if len(target_to_query[id2])<10:
        continue
    if id2 not in target_data:
        continue

    ref_data = (id2, target_data[id2]) 
    ref_len = len(target_data[id2])
    base_count.clear()
    base_count.update( target_data[id2] )
    if 1.0*base_count.most_common(1)[0][1]/ref_len > 0.8:  # don't do preassmbly if a read is of >80% of the same base
        continue
    read_data = []
    
    query_alignment = target_to_query[id2]
    query_alignment.sort() # get better alignment
    total_bases = 0
    max_cov_bases = max_cov * ref_len * 1.2
    #min_cov_bases = min_cov * ref_len * 3
    
    for rank_score, l in query_alignment:
        rank, score = rank_score
        #score = rank_score
        l = l.split()
        id1 = l[0]
        #if -noSplitSubread not used, we will need the following line    
        #id1 = id1.split("/")[0]
        q_s = int(l[5]) + trim_align
        q_e = int(l[6]) - trim_align
        strand = int(l[8])
        t_s = int(l[9])
        t_e = int(l[10])
        t_l = int(l[11])
        #if strand == 1:
        #    t_s, t_e = t_l - t_e, t_l - t_s
        #    t_s += trim_align
        #    t_e -= trim_align

        if q_e - q_s < 400:
            continue
        total_bases += q_e - q_s
        if total_bases > max_cov_bases:
            break
        q_seq = query_data[id1][q_s:q_e]
        read_data.append( ( "%s/0/%d_%d" % (id1, q_s, q_e), q_s, q_e, q_seq, strand, t_s, t_e) )

    if len(read_data) > 5:
        r_count += 1
        t_id, t_seq = ref_data 
        t_len = len(t_seq)
        print t_id, t_seq
        for r in read_data:
            q_id, q_s, q_e, q_seq, strand, t_s, t_e = r
            if strand == 1:
                q_seq = rc(q_seq)
            print q_id, q_seq
        #if r_count > 600:
        #    break
        print "+ +"
print "- -"

#output_dir,dumb = os.path.split( os.path.abspath( output_file ) )
#output_log = open ( os.path.join( output_dir, "j%02d.log" % group_id ), "w" )




