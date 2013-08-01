from ctypes import *
import sys
from multiprocessing import Pool

falcon = CDLL("./falcon.so")

falcon.generate_consensus.argtypes = [POINTER(c_char_p), c_uint ]
falcon.generate_consensus.restype = POINTER(c_char)
falcon.free_consensus.argtypes = [ c_char_p ]

def get_consensus( c_input ):
    seqs, seed_id = c_input
    seqs_ptr = (c_char_p * len(seqs))()
    seqs_ptr[:] = seqs
    consensus_ptr = falcon.generate_consensus( seqs_ptr, len(seqs) )
    consensus = string_at(consensus_ptr)[:]
    falcon.free_consensus( consensus_ptr )
    del seqs_ptr
    return consensus, seed_id

exe_pool = Pool(24)


def get_seq_data():
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
                    yield (seqs, seed_id) 
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
