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
import os


rcmap = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))

if __name__ == "__main__":
    import argparse
    import re
    from pbcore.io import FastaReader
    
    tiling_path = {}
    with open("primary_tigs_tiling_path_c") as f:
        for l in f:
            l = l.strip().split()
            tiling_path.setdefault( l[0], [])

            offset = int(l[1])
            node_id = l[2].split(":")
            s = int(l[3])
            e = int(l[4])

            tiling_path[ l[0] ].append( (offset, node_id[0], node_id[1], s, e) )

    f = FastaReader("preads.fa")
    seq_db = {}
    for r in f:
         seq_db[r.name] = r.sequence

    f = FastaReader("primary_tigs_c.fa")
    ptigs_db = {}
    for r in f:
         ptigs_db[r.name] = r.sequence

    for ptig_id in ptigs_db:
        pread_data = {}
        offsets = []
        seqs = []
        p_tig = ptigs_db[ptig_id]
        print ptig_id, 0, p_tig
        for offset, s_id, end, s, e in tiling_path[ptig_id]:
            seq = seq_db[s_id]
            if end == "B":
                s, e = e, s
                offset = offset - len(seq) 
                seq = "".join([rcmap[c] for c in seq[::-1]])
            else:
                offset = offset - len(seq)
            print s_id, offset, seq
        
        print "+ + +"
    print "- - -"

