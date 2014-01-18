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

from ctypes import *
import sys
from multiprocessing import Pool
import os
import falcon_kit
module_path = falcon_kit.__path__[0]

falcon = CDLL(os.path.join(module_path, "falcon.so"))

falcon.generate_consensus.argtypes = [POINTER(c_char_p), c_uint ]
falcon.generate_consensus.restype = POINTER(c_char)
falcon.free_consensus.argtypes = [ c_char_p ]

def get_consensus( c_input ):
    seqs, seed_id, min_cov, K = c_input
    seqs_ptr = (c_char_p * len(seqs))()
    seqs_ptr[:] = seqs
    consensus_ptr = falcon.generate_consensus( seqs_ptr, len(seqs), min_cov, K )
    consensus = string_at(consensus_ptr)[:]
    falcon.free_consensus( consensus_ptr )
    del seqs_ptr
    return consensus, seed_id

exe_pool = Pool(24)


def get_seq_data(min_cov = 8, K = 8):
    seqs = []
    seed_id = None
    seqs_data = []
    with sys.stdin as f:
        for l in f:
            l = l.strip().split()
            if l[0] not in ("+", "-"):
                if len(seqs) == 0:
                    seqs.append(l[1]) #the "seed"
                    seed_id = l[0]
                seqs.append(l[1])
            elif l[0] == "+":
                if len(seqs) > 10:
                    yield (seqs, seed_id, min_cov, K) 
                #seqs_data.append( (seqs, seed_id) ) 
                seqs = []
                seed_id = None
            elif l[0] == "-":
                #yield (seqs, seed_id)
                #seqs_data.append( (seqs, seed_id) )
                break

for res in exe_pool.imap(get_consensus, get_seq_data()):
    cns, seed_id = res
    if len(cns) > 500:
        print ">"+seed_id
        print cns
