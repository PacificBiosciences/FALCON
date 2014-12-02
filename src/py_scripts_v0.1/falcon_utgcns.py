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
"""
consensus_data * generate_utg_consensus( char ** input_seq, 
                           seq_coor_t *offset,
                           unsigned int n_seq, 
                           unsigned min_cov, 
                           unsigned K,
                           double min_idt) {
"""
falcon.generate_utg_consensus.argtypes = [ POINTER(c_char_p), POINTER(falcon_kit.seq_coor_t), c_uint, c_uint, c_uint, c_double ]
falcon.generate_utg_consensus.restype = POINTER(falcon_kit.ConsensusData)
falcon.free_consensus_data.argtypes = [ POINTER(falcon_kit.ConsensusData) ]

rcmap = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))

def get_consensus(c_input):
    t_id, seqs, offsets, config = c_input 
    K = config[0]
    seqs_ptr = (c_char_p * len(seqs))()
    seqs_ptr[:] = seqs
    offset_ptr = (c_long * len(seqs))( *offsets )
    consensus_data_ptr = falcon.generate_utg_consensus( seqs_ptr, offset_ptr, len(seqs), 0, K, 0.)
    consensus = string_at(consensus_data_ptr[0].sequence)[:]
    del seqs_ptr
    del offset_ptr
    falcon.free_consensus_data( consensus_data_ptr )
    return consensus, t_id

def echo(c_input):

    t_id, seqs, offsets, config = c_input 

    return len(seqs), "test"

def get_seq_data(config):
    seqs = []
    offsets = []
    seed_id = None
    with sys.stdin as f:
        for l in f:
            l = l.strip().split()
            if len(l) != 3:
                continue
            if l[0] not in ("+", "-"):
                if len(seqs) == 0:
                    seqs.append(l[2]) #the "seed"
                    offsets.append( int(l[1]) )
                    seed_id = l[0]
                else:
                    seqs.append(l[2])
                    offsets.append( int(l[1]) )
            elif l[0] == "+":
                yield (seed_id, seqs, offsets, config) 
                seqs = []
                offsets = []
                seed_id = None
            elif l[0] == "-":
                break

if __name__ == "__main__":
    import argparse
    import re
    parser = argparse.ArgumentParser(description='a simple multi-processor consensus sequence generator')
    parser.add_argument('--n_core', type=int, default=4,
                        help='number of processes used for generating consensus')
    args = parser.parse_args()
    exe_pool = Pool(args.n_core)
    K = 8
    config = (K, )
    for res in exe_pool.imap(get_consensus, get_seq_data(config)):  
    #for res in exe_pool.imap(echo, get_seq_data(config)):  
    #for res in map(echo, get_seq_data(config)):  
    #for res in map(get_consensus, get_seq_data(config)):  
        cns, t_id = res
        print ">"+t_id+"|tigcns"
        print cns

