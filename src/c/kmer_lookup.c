/*
 * =====================================================================================
 *
 *       Filename:  kmer_count.c
 *
 *    Description:  
 *
 *        Version:  0.1
 *        Created:  07/20/2013 17:00:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jason Chin, 
 *        Company:  
 *
 * =====================================================================================

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
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "common.h"


const unsigned int KMERMATCHINC = 10000;


kmer_lookup * allocate_kmer_lookup ( seq_coor_t size ) {
    kmer_lookup * kl;
    seq_coor_t i;

    //printf("%lu is allocated for kmer lookup\n", size);
    kl = (kmer_lookup *)  malloc( size * sizeof(kmer_lookup) );
    init_kmer_lookup( kl, size);
    return kl;
}

void init_kmer_lookup ( kmer_lookup * kl,  seq_coor_t size ) {
    seq_coor_t i;
    //printf("%lu is allocated for kmer lookup\n", size);
    for (i=0; i<size; i++) {
        kl[i].start = LONG_MAX;
        kl[i].last = LONG_MAX;
        kl[i].count = 0;
    }
}


void free_kmer_lookup( kmer_lookup *  ptr) {
    free(ptr);
}

seq_array allocate_seq(seq_coor_t size) {
    seq_array sa;
    sa  = (seq_array) malloc( size * sizeof(base) ); 
    init_seq_array( sa, size);
    return sa;
}

void init_seq_array( seq_array sa, seq_coor_t size) {
    seq_coor_t i;
    for (i=0; i++; i<size) {
        sa[i] = 0xff;
    }
}

void free_seq_array( seq_array sa) {
    free(sa);
}

seq_addr_array allocate_seq_addr(seq_coor_t size) {
    return (seq_addr_array) calloc( size, sizeof(seq_addr));
}

void free_seq_addr_array(seq_addr_array sda) {
    free(sda);
}

seq_coor_t get_kmer_bitvector(seq_array sa, unsigned int K) {
    unsigned int i;
    seq_coor_t kmer_bv = 0;
    seq_coor_t kmer_mask;

    kmer_mask = 0;
    for (i = 0; i < K; i++) {
        kmer_mask <<= 2;
        kmer_mask |= 0x00000003;
    }

    for (i = 0; i < K; i++) {
        kmer_bv <<= 2;
        kmer_bv |= (unsigned int) sa[i];
    }

    return kmer_bv;
}

void add_sequence ( seq_coor_t start, 
                    unsigned int K, 
                    char * seq, 
                    seq_coor_t seq_len,
                    seq_addr_array sda, 
                    seq_array sa, 
                    kmer_lookup * lk ) {

    seq_coor_t i;
    seq_coor_t kmer_bv;
    seq_coor_t kmer_mask;

    kmer_mask = 0;
    for (i = 0; i < K; i++) {
        kmer_mask <<= 2;
        kmer_mask |= 0x00000003;
    }

    for (i = 0; i < seq_len; i++) {
        switch ( seq[i] ) {
            case 'A':
                sa[ start + i ] = 0;
                break;
            case 'C':
                sa[ start + i ] = 1;
                break;
            case 'G':
                sa[ start + i ] = 2;
                break;
            case 'T':
                sa[ start + i ] = 3;
        }
    }
    kmer_bv = get_kmer_bitvector( sa + start, K);
    for (i = 0; i < seq_len - K;  i++) {
        //printf("%lu %lu\n", i, kmer_bv);
        //printf("lk before init: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
        if (lk[kmer_bv].start == LONG_MAX) {
            lk[kmer_bv].start = start + i;
            lk[kmer_bv].last = start + i;
            lk[kmer_bv].count += 1;
            //printf("lk init: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
        } else {
            sda[ lk[kmer_bv].last ] = start + i;
            lk[kmer_bv].count += 1;
            lk[kmer_bv].last = start + i;
            //printf("lk change: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
        }
        kmer_bv <<= 2;
        kmer_bv |= sa[ start + i + K];
        kmer_bv &= kmer_mask;
    }
}


void mask_k_mer(seq_coor_t size, kmer_lookup * kl, seq_coor_t threshold) {
    seq_coor_t i;
    for (i=0; i<size; i++) {
        if (kl[i].count > threshold) {
            kl[i].start = LONG_MAX;
            kl[i].last = LONG_MAX;
            //kl[i].count = 0;
        }
    }
}


kmer_match * find_kmer_pos_for_seq( char * seq, seq_coor_t seq_len, unsigned int K,
                    seq_addr_array sda, 
                    kmer_lookup * lk) {
    seq_coor_t i;
    seq_coor_t kmer_bv;
    seq_coor_t kmer_mask;
    seq_coor_t kmer_pos;
    seq_coor_t next_kmer_pos;
    unsigned int half_K;
    seq_coor_t kmer_match_rtn_allocation_size = KMERMATCHINC;
    kmer_match * kmer_match_rtn;
    base * sa;

    kmer_match_rtn = (kmer_match *) malloc( sizeof(kmer_match) );
    kmer_match_rtn->count = 0;
    kmer_match_rtn->query_pos = (seq_coor_t *) calloc( kmer_match_rtn_allocation_size, sizeof( seq_coor_t ) );
    kmer_match_rtn->target_pos = (seq_coor_t *) calloc( kmer_match_rtn_allocation_size, sizeof( seq_coor_t ) );

    sa = calloc( seq_len, sizeof(base) );

    kmer_mask = 0;
    for (i = 0; i < K; i++) {
        kmer_mask <<= 2;
        kmer_mask |= 0x00000003;
    }

    for (i = 0; i < seq_len; i++) {
        switch ( seq[i] ) {
            case 'A':
                sa[ i ] = 0;
                break;
            case 'C':
                sa[ i ] = 1;
                break;
            case 'G':
                sa[ i ] = 2;
                break;
            case 'T':
                sa[ i ] = 3;
        }
    }


    kmer_bv = get_kmer_bitvector(sa, K);
    half_K = K >> 1;
    for (i = 0; i < seq_len - K;  i += half_K) {
        kmer_bv = get_kmer_bitvector(sa + i, K);
        if (lk[kmer_bv].start == LONG_MAX) {  //for high count k-mers
            continue;
        }
        kmer_pos = lk[ kmer_bv ].start;
        next_kmer_pos = sda[ kmer_pos ];
        kmer_match_rtn->query_pos[ kmer_match_rtn->count ] = i;
        kmer_match_rtn->target_pos[ kmer_match_rtn->count ] = kmer_pos;
        kmer_match_rtn->count += 1;
        if (kmer_match_rtn->count > kmer_match_rtn_allocation_size - 1000) {
            kmer_match_rtn_allocation_size += KMERMATCHINC;
            kmer_match_rtn->query_pos = (seq_coor_t *) realloc( kmer_match_rtn->query_pos, 
                                                                   kmer_match_rtn_allocation_size  * sizeof(seq_coor_t) );
            kmer_match_rtn->target_pos = (seq_coor_t *) realloc( kmer_match_rtn->target_pos, 
                                                                    kmer_match_rtn_allocation_size  * sizeof(seq_coor_t) );
        }
        while ( next_kmer_pos > kmer_pos ){
            kmer_pos = next_kmer_pos;
            next_kmer_pos = sda[ kmer_pos ];
            kmer_match_rtn->query_pos[ kmer_match_rtn->count ] = i;
            kmer_match_rtn->target_pos[ kmer_match_rtn->count ] = kmer_pos;
            kmer_match_rtn->count += 1;
            if (kmer_match_rtn->count > kmer_match_rtn_allocation_size - 1000) {
                kmer_match_rtn_allocation_size += KMERMATCHINC;
                kmer_match_rtn->query_pos = (seq_coor_t *) realloc( kmer_match_rtn->query_pos, 
                                                                       kmer_match_rtn_allocation_size  * sizeof(seq_coor_t) );
                kmer_match_rtn->target_pos = (seq_coor_t *) realloc( kmer_match_rtn->target_pos, 
                                                                        kmer_match_rtn_allocation_size  * sizeof(seq_coor_t) );
            }
        }
    }
    free(sa);
    return kmer_match_rtn;
}

void free_kmer_match( kmer_match * ptr) {
    free(ptr->query_pos);
    free(ptr->target_pos);
    free(ptr);
}

aln_range* find_best_aln_range(kmer_match * km_ptr, 
                              seq_coor_t K, 
                              seq_coor_t bin_size, 
                              seq_coor_t count_th) {
    seq_coor_t i;
    seq_coor_t j;
    seq_coor_t q_min, q_max, t_min, t_max;
    seq_coor_t * d_count;
    seq_coor_t * q_coor;
    seq_coor_t * t_coor;
    aln_range * arange;

    long int d, d_min, d_max;
    long int cur_score;
    long int max_score;
    long int max_k_mer_count;
    long int max_k_mer_bin;
    seq_coor_t cur_start;
    seq_coor_t cur_pos;
    seq_coor_t max_start;
    seq_coor_t max_end;
    seq_coor_t kmer_dist;

    arange = calloc(1 , sizeof(aln_range));

    q_min = LONG_MAX;
    q_max = 0;
    t_min = LONG_MAX;
    t_max = 0;

    d_min = LONG_MAX;
    d_max = LONG_MIN;

    for (i = 0; i <  km_ptr->count; i++ ) {
        if ( km_ptr -> query_pos[i] < q_min) {
            q_min =  km_ptr->query_pos[i];
        }
        if ( km_ptr -> query_pos[i] > q_max) {
            q_max =  km_ptr->query_pos[i];
        }
        if ( km_ptr -> target_pos[i] < t_min) {
            t_min =  km_ptr->target_pos[i];
        }
        if ( km_ptr -> query_pos[i] > t_max) {
            t_max =  km_ptr->target_pos[i];
        }
        d = (long int) km_ptr->query_pos[i] - (long int) km_ptr->target_pos[i];
        if ( d < d_min ) {
            d_min = d;
        }
        if ( d > d_max ) {
            d_max = d;
        }
    }

    //printf("%lu %ld %ld\n" , km_ptr->count, d_min, d_max);
    d_count = calloc( (d_max - d_min)/bin_size + 1, sizeof(seq_coor_t) );
    q_coor = calloc( km_ptr->count, sizeof(seq_coor_t) );
    t_coor = calloc( km_ptr->count, sizeof(seq_coor_t) );

    for (i = 0; i <  km_ptr->count; i++ ) {
        d = (long int) (km_ptr->query_pos[i]) - (long int) (km_ptr->target_pos[i]);
        d_count[ (d - d_min)/ (long int) bin_size ] += 1;
        q_coor[i] = LONG_MAX;
        t_coor[i] = LONG_MAX;
    }

    j = 0;
    max_k_mer_count = 0;
    max_k_mer_bin = LONG_MAX;
    for (i = 0; i <  km_ptr->count; i++ ) {
        d = (long int) (km_ptr->query_pos[i]) - (long int) (km_ptr->target_pos[i]);
        if ( d_count[ (d - d_min)/ (long int) bin_size ] > max_k_mer_count) {
            max_k_mer_count =  d_count[ (d - d_min)/ (long int) bin_size ];
            max_k_mer_bin = (d - d_min)/ (long int) bin_size;
        }
    }
    //printf("k_mer: %lu %lu\n" , max_k_mer_count, max_k_mer_bin);
    
    if ( max_k_mer_bin != LONG_MAX && max_k_mer_count > count_th ) {
        for (i = 0; i <  km_ptr->count; i++ ) {
            d = (long int) (km_ptr->query_pos[i]) - (long int) (km_ptr->target_pos[i]);
            if ( abs( ( (d - d_min)/ (long int) bin_size ) - max_k_mer_bin ) > 10 ) {
                continue;
            }
            if (d_count[ (d - d_min)/ (long int) bin_size ] > count_th) {
                q_coor[j] = km_ptr->query_pos[i];  
                t_coor[j] = km_ptr->target_pos[i];
                //printf("d_count: %lu %lu\n" ,i, d_count[(d - d_min)/ (long int) bin_size]);
                //printf("coor: %lu %lu\n" , q_coor[j], t_coor[j]);
                j ++;
            }
        }
    }

    if (j > 1) {
        arange->s1 = q_coor[0];
        arange->e1 = q_coor[0];
        arange->s2 = t_coor[0];
        arange->e2 = t_coor[0];

        max_score = 0;
        cur_score = 0;
        cur_start = 0;

        for (i = 1; i < j; i++) {
            cur_score += 32 - (q_coor[i] - q_coor[i-1]);
            //printf("deltaD, %lu %ld\n", q_coor[i] - q_coor[i-1], cur_score);
            if (cur_score < 0) {
                cur_score = 0;
                cur_start = i;
            } else if (cur_score > max_score) {
                arange->s1 = q_coor[cur_start];
                arange->s2 = t_coor[cur_start];
                arange->e1 = q_coor[i];
                arange->e2 = t_coor[i];
                max_score = cur_score;
                //printf("%lu %lu %lu %lu\n", arange.s1, arange.e1, arange.s2, arange.e2);
            }
        }

    } else {
        arange->s1 = 0;
        arange->e1 = 0;
        arange->s2 = 0;
        arange->e2 = 0;
    }

    // printf("free\n");

    free(d_count);
    free(q_coor);
    free(t_coor);
    return arange;
}

void free_aln_range( aln_range * arange) {
    free(arange);
}
