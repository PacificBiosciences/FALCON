from ctypes import *
from pbcore.io import FastaReader
import sys
from multiprocessing import Pool

falcon = CDLL("./falcon.so")

falcon.generate_consensus.argtypes = [POINTER(c_char_p), c_uint ]
falcon.generate_consensus.restype = POINTER(c_char)
falcon.free_consensus.argtypes = [ c_char_p ]
RCMAP = dict(zip("ACGTacgtNn-","TGCAtgcaNn-"))

def get_consensus( c_input ):
    seqs, seed_id, min_cov, K = c_input
    seqs_ptr = (c_char_p * len(seqs))()
    seqs_ptr[:] = seqs
    consensus_ptr = falcon.generate_consensus( seqs_ptr, len(seqs), min_cov, K )
    consensus = string_at(consensus_ptr)[:]
    falcon.free_consensus( consensus_ptr )
    del seqs_ptr
    return consensus, seed_id


def get_seq_data(min_cov = 0, K = 14):
    reads = []
    for r in FastaReader(sys.argv[1]):
        reads.append(r.sequence)
        reads.append("".join([RCMAP[c] for c in r.sequence[::-1]]))


    for r in FastaReader(sys.argv[2]):
        ref_seq = r.sequence
        ref_name = r.name
        seqs = [r.sequence]
        seqs.extend(reads)
        yield (seqs, ref_name, min_cov, K) 
        

exe_pool = Pool(4)
for res in exe_pool.imap(get_consensus, get_seq_data()):
    cns, seed_id = res
    if len(cns) > 500:
        print ">"+seed_id
        print cns
