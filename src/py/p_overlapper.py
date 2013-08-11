from falcon_kit import * 
from pbcore.io import FastaReader
import numpy as np
import collections
import sys
import multiprocessing as mp
from multiprocessing import sharedctypes
from ctypes import *

global sa_ptr, sda_ptr, lk_ptr
global q_seqs, seqs
RC_MAP = dict( zip("ACGTacgtNn-", "TGCAtgcaNn-") )

def get_ovelap_alignment(seq1, seq0):

    K = 8
    lk_ptr = kup.allocate_kmer_lookup( 1 << (K * 2) )
    sa_ptr = kup.allocate_seq( len(seq0) )
    sda_ptr = kup.allocate_seq_addr( len(seq0) )
    kup.add_sequence( 0, K, seq0, len(seq0), sda_ptr, sa_ptr, lk_ptr)

    kmer_match_ptr = kup.find_kmer_pos_for_seq(seq1, len(seq1), K, sda_ptr, lk_ptr)
    kmer_match = kmer_match_ptr[0]
    aln_range_ptr = kup.find_best_aln_range(kmer_match_ptr, K, K*5, 50)
    #x,y = zip( * [ (kmer_match.query_pos[i], kmer_match.target_pos[i]) for i in range(kmer_match.count )] )
    aln_range = aln_range_ptr[0]
    kup.free_kmer_match(kmer_match_ptr)
    s1, e1, s0, e0 = aln_range.s1, aln_range.e1, aln_range.s2, aln_range.e2  
    e1 += K + K/2
    e0 += K + K/2
    kup.free_aln_range(aln_range)
    len_1 = len(seq1)
    len_0 = len(seq0)
    if e1 > len_1: 
        e1 = len_1
    if e0 > len_0:
        e0 = len_0
    do_aln = False
    contain_status = "none" 
    #print s0, e0, s1, e1 
    if e1 - s1 > 500:
        if s0 < s1 and s0 > 24:
            do_aln = False
        elif s1 <= s0 and s1 > 24:
            do_aln = False
        elif s1 < 24 and len_1 - e1 < 24:
            do_aln = False
            contain_status = "contains"
            #print "X1"
        elif s0 < 24 and len_0 - e0 < 24:
            do_aln = False
            contain_status = "contained"
            #print "X2"
        else:
            do_aln = True
            if s0 < s1:
                s1 -= s0 #assert s1 > 0
                s0 = 0
                e1 = len_1
                #if len_1 - s1 >= len_0:
                #    do_aln = False
                #    contain_status = "contains"
                #    print "X3", s0, e0, len_0, s1, e1, len_1

                
            elif s1 <= s0:
                s0 -= s1 #assert s1 > 0
                s1 = 0
                e0 = len_0
                print s0, e0, s1, e1
                #if len_0 - s0 >= len_1:
                #    do_aln = False
                #    contain_status = "contained"
                #    print "X4"
        #if abs( (e1 - s1) - (e0 - s0 ) ) > 200:  #avoid overlap alignment for big indels
        #    do_aln = False

        if do_aln:
            alignment = DWA.align(seq1[s1:e1], e1-s1,
                                  seq0[s0:e0], e0-s0,
                                  500, 0)
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
            assert aln_q_e- aln_q_s <= alignment[0].aln_str_size or aln_t_e- aln_t_s <= alignment[0].aln_str_size
            #print aln_str1
            #print aln_str0
            if aln_size > 500: 
                contain_status = "overlap"            
            DWA.free_alignment(alignment)
        
    kup.free_seq_addr_array(sda_ptr)
    kup.free_seq_array(sa_ptr)
    kup.free_kmer_lookup(lk_ptr)

    if e1 - s1 > 500 and do_aln and aln_size > 500:
        #return s1, s1+aln_q_e-aln_q_s, s2, s2+aln_t_e-aln_t_s, aln_size, aln_dist, x, y
        return s1, s1+aln_q_e-aln_q_s, s0, s0+aln_t_e-aln_t_s, aln_size, aln_dist, contain_status
    else:
        return 0, 0, 0, 0, 0, 0, contain_status 

def get_candidate_hits(q_name):
    global sa_ptr, sda_ptr, lk_ptr
    global q_seqs
    K = 14
    q_seq = q_seqs[q_name]

    rtn = []

    c_sda_ptr = cast(sda_ptr, POINTER(seq_coor_t))
    c_sa_ptr = cast(sa_ptr, POINTER(base_t))
    c_lk_ptr = cast(lk_ptr, POINTER(KmerLookup))

    kmer_match_ptr = kup.find_kmer_pos_for_seq(q_seq, len(q_seq), K, c_sda_ptr, c_lk_ptr)
    kmer_match = kmer_match_ptr[0]
    count = kmer_match.count
    hit_index = np.array(kmer_match.target_pos[0:count])/500
    kup.free_kmer_match(kmer_match_ptr)
    
    c = collections.Counter(hit_index)
    s = [c[0] for c in c.items() if c[1] >50]
    #s.sort()
    targets = set()
    for p in s:
        hit_id = seqs[p/2][0]
        if hit_id in targets or hit_id == q_name:
            continue
        targets.add(hit_id)
        seq1, seq0 = q_seq, q_seqs[hit_id]
        aln_data = get_ovelap_alignment(seq1, seq0)
        #rtn = get_alignment(seq1, seq0)
        if rtn != None:
            
            s1, e1, s2, e2, aln_size, aln_dist, c_status = aln_data
            #print >>f, name, 0, s1, e1, len(seq1), hit_id, 0, s2, e2, len(seq0),  aln_size, aln_dist
            rtn.append( ( hit_id, q_name, aln_dist - aln_size, "%0.2f" % (100 - 100.0*aln_dist/(aln_size+1)), 
                          0, s2, e2, len(seq0), 
                          0, s1, e1, len(seq1), c_status ) )

    r_q_seq = "".join([RC_MAP[c] for c in q_seq[::-1]])
    
    kmer_match_ptr = kup.find_kmer_pos_for_seq(r_q_seq, len(r_q_seq), K, c_sda_ptr, c_lk_ptr)
    kmer_match = kmer_match_ptr[0]
    count = kmer_match.count
    hit_index = np.array(kmer_match.target_pos[0:count])/500
    kup.free_kmer_match(kmer_match_ptr)
    
    c = collections.Counter(hit_index)
    s = [c[0] for c in c.items() if c[1] >50]
    #s.sort()
    targets = set()
    for p in s:
        hit_id = seqs[p/2][0]
        if hit_id in targets or hit_id == q_name:
            continue
        targets.add(hit_id)
        seq1, seq0 = r_q_seq, q_seqs[hit_id]
        aln_data = get_ovelap_alignment(seq1, seq0)
        #rtn = get_alignment(seq1, seq0)
        if rtn != None:
            s1, e1, s2, e2, aln_size, aln_dist, c_status = aln_data
            #print >>f, name, 1, s1, e1, len(seq1), hit_id, 0, s2, e2, len(seq0),  aln_size, aln_dist
            rtn.append( ( hit_id, q_name, aln_dist - aln_size, "%0.2f" % (100 - 100.0*aln_dist/(aln_size+1)), 
                          0, s2, e2, len(seq0), 
                          1, len(seq1) - e1, len(seq1)- s1, len(seq1), c_status ) )

    return rtn

def build_look_up(seqs, K):
    global sa_ptr, sda_ptr, lk_ptr

    total_index_base = len(seqs) * 1000
    sa_ptr = sharedctypes.RawArray(base_t, total_index_base)
    c_sa_ptr = cast(sa_ptr, POINTER(base_t))
    kup.init_seq_array(c_sa_ptr, total_index_base)

    sda_ptr = sharedctypes.RawArray(seq_coor_t, total_index_base)
    c_sda_ptr = cast(sda_ptr, POINTER(seq_coor_t))

    lk_ptr = sharedctypes.RawArray(KmerLookup, 1 << (K*2))
    c_lk_ptr = cast(lk_ptr, POINTER(KmerLookup))
    kup.init_kmer_lookup(c_lk_ptr, 1 << (K*2))

    start = 0
    for r_name, prefix, suffix in seqs:
        kup.add_sequence( start, K, prefix, 500, c_sda_ptr, c_sa_ptr, c_lk_ptr)
        start += 500
        kup.add_sequence( start, K, suffix, 500, c_sda_ptr, c_sa_ptr, c_lk_ptr)
        start += 500

    kup.mask_k_mer(1 << (K * 2), c_lk_ptr, 256)
    
    #return sda_ptr, sa_ptr, lk_ptr


def lookup_data_iterator( q_seqs ):
    for q_name, q_seq in q_seqs.items():
        yield q_name


if __name__ == "__main__":
    seqs = []
    q_seqs = {}
    f = FastaReader(sys.argv[1]) # take one commnad line argument of the input fasta file name

    for r in f:
        if len(r.sequence) < 6000:
            continue
        seq = r.sequence.upper()
        seqs.append( (r.name, seq[:500], seq[-500:] ) )
        q_seqs[r.name] = seq


    total_index_base = len(seqs) * 1000
    K = 14
    build_look_up(seqs, K)
    pool = mp.Pool(8)

    for r in pool.imap(get_candidate_hits, lookup_data_iterator( q_seqs)):
        for h in r:
            print " ".join([str(x) for x in h]) 

