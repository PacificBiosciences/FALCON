from ctypes import *


seq_coor_t = c_long

class KmerLookup(Structure):
    _fields_ = [("start", seq_coor_t),
                ("last", seq_coor_t),
                ("count", seq_coor_t)]

class KmerMatch(Structure):
    _fields_ = [ ("count", seq_coor_t),
                ("query_pos", POINTER(seq_coor_t)),
                ("target_pos", POINTER(seq_coor_t)) ]

class AlnRange(Structure):
    _fields_ = [ ("s1", seq_coor_t),
                 ("e1", seq_coor_t),
                 ("s2", seq_coor_t),
                 ("e2", seq_coor_t) ]

kup = CDLL("./kmer_lookup.so")

kup.allocate_kmer_lookup.argtypes =  [seq_coor_t] 
kup.allocate_kmer_lookup.restype = POINTER(KmerLookup)
kup.free_kmer_lookup.argtypes = [POINTER(KmerLookup)]

kup.allocate_seq.argtypes = [seq_coor_t]
kup.allocate_seq.restype = POINTER(c_uint8)
kup.free_seq_array.argtypes = [POINTER(c_uint8)]

kup.allocate_seq_addr.argtypes = [seq_coor_t]
kup.allocate_seq_addr.restype = POINTER(seq_coor_t)
kup.free_seq_addr_array.argtypes = [POINTER(seq_coor_t)]

kup.add_sequence.argtypes = [ seq_coor_t, c_uint, POINTER(c_char), seq_coor_t, POINTER(seq_coor_t), POINTER(c_uint8), POINTER(KmerLookup) ]
kup.mask_k_mer.argtypes =[ c_long, POINTER(KmerLookup), c_long ]
kup.find_kmer_pos_for_seq.argtypes = [ POINTER(c_char), seq_coor_t, c_uint, POINTER(seq_coor_t), POINTER(KmerLookup)]
kup.find_kmer_pos_for_seq.restype = POINTER(KmerMatch)
kup.free_kmer_match.argtypes = [ POINTER(KmerMatch) ]


kup.find_best_aln_range.argstypes = [POINTER(KmerMatch), seq_coor_t, seq_coor_t, seq_coor_t]
kup.find_best_aln_range.restype = AlnRange


class Alignment(Structure):
    """
    typedef struct {    
        seq_coor_t aln_str_size ;
        seq_coor_t dist ;
        seq_coor_t aln_q_s;
        seq_coor_t aln_q_e;
        seq_coor_t aln_t_s;
        seq_coor_t aln_t_e;
        char * q_aln_str;
        char * t_aln_str;
    } alignment;
    """
    _fields_ = [ ("aln_str_size", seq_coor_t),
                 ("dist", seq_coor_t),
                 ("aln_q_s", seq_coor_t),
                 ("aln_q_e", seq_coor_t),
                 ("aln_t_s", seq_coor_t),
                 ("aln_t_e", seq_coor_t),
                 ("q_aln_str", c_char_p),
                 ("t_aln_str", c_char_p)]


DWA = CDLL("./DW_align.so")
DWA.align.argtypes = [ POINTER(c_char), c_long, POINTER(c_char), c_long, c_long, c_int ] 
DWA.align.restype = POINTER(Alignment)
DWA.free_alignment.argtypes = [POINTER(Alignment)]



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



def get_alignment(seq1, seq0):
    K = 8
    lk_ptr = kup.allocate_kmer_lookup( 1 << (K * 2) )
    sa_ptr = kup.allocate_seq( len(seq0) )
    sda_ptr = kup.allocate_seq_addr( len(seq0) )
    kup.add_sequence( 0, K, seq0, len(seq0), sda_ptr, sa_ptr, lk_ptr)

    kmer_match_ptr = kup.find_kmer_pos_for_seq(seq1, len(seq1), K, sda_ptr, lk_ptr)
    kmer_match = kmer_match_ptr[0]
    aln_range = kup.find_best_aln_range(kmer_match_ptr, K, K*10, 50)
    #x,y = zip( * [ (kmer_match.query_pos[i], kmer_match.target_pos[i]) for i in range(kmer_match.count )] )
    kup.free_kmer_match(kmer_match_ptr)
    s1, e1, s2, e2 = aln_range.s1, aln_range.e1, aln_range.s2, aln_range.e2
    if e1 - s1 > 500:
        #s1 = 0 if s1 < 14 else s1 - 14
        #s2 = 0 if s2 < 14 else s2 - 14
        e1 = len(seq1) if e1 >= len(seq1)-2*K else e1 + K*2
        e2 = len(seq0) if e2 >= len(seq0)-2*K else e2 + K*2
        
        alignment = DWA.align(seq1[s1:e1], e1-s1,
                              seq0[s2:e2], e2-s2,
                              100,
                              0)
        #print seq1[s1:e1]
        #print seq0[s2:e2]
        #if alignment[0].aln_str_size > 500:

        #aln_str1 = alignment[0].q_aln_str
        #aln_str0 = alignment[0].t_aln_str
        aln_size = alignment[0].aln_str_size
        aln_dist = alignment[0].dist
        aln_q_s = alignment[0].aln_q_s
        aln_q_e = alignment[0].aln_q_e
        aln_t_s = alignment[0].aln_t_s
        aln_t_e = alignment[0].aln_t_e
        
        #print "X,",alignment[0].aln_q_s, alignment[0].aln_q_e
        #print "Y,",alignment[0].aln_t_s, alignment[0].aln_t_e
        
        #print aln_str1
        #print aln_str0
    
        DWA.free_alignment(alignment)

    kup.free_seq_addr_array(sda_ptr)
    kup.free_seq_array(sa_ptr)
    kup.free_kmer_lookup(lk_ptr)
    if e1 - s1 > 500 and aln_size > 500:
        return s1, s1+aln_q_e-aln_q_s, s2, s2+aln_t_e-aln_t_s, aln_size, aln_dist
    else:
        return None
