/*
 * =====================================================================================
 *
 *       Filename:  DW_banded.c
 *
 *    Description:  A banded version for the O(ND) greedy sequence alignment algorithm 
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
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <stdbool.h>
#include "common.h"


alignment * align(char * query_seq, seq_coor_t q_len,
                  char * target_seq, seq_coor_t t_len,
                  seq_coor_t band_tolerance,
                  int get_aln_str) {
    seq_coor_t * V;
    seq_coor_t * U;  // array of matched bases for each "k"
    seq_coor_t k_offset;
    seq_coor_t d;
    seq_coor_t k, k2;
    seq_coor_t best_m;  // the best "matches" for each d
    seq_coor_t min_k, new_min_k;
    seq_coor_t max_k, new_max_k;
    seq_coor_t pre_k;
    seq_coor_t x, y;
    seq_coor_t cd;
    seq_coor_t ck;
    seq_coor_t cx, cy, nx, ny;
    seq_coor_t max_d;
    seq_coor_t band_size;

    d_path_data * d_path;
    path_point * aln_path;
    seq_coor_t aln_path_idx;
    alignment * align_rtn;
    seq_coor_t aln_pos;
    seq_coor_t i;
    bool aligned = false;

    //printf("debug: %ld %ld\n", q_len, t_len);
    //printf("%s\n", query_seq);
   
    max_d = q_len + t_len;

    //band_size = max_d > 5000 ? 5000 : max_d;
    band_size = (max_d >> 2) > 600 ? 600 : (max_d >> 2);

    V = calloc( (band_size + 1) * 2 + 1, sizeof(seq_coor_t) );
    U = calloc( (band_size + 1) * 2 + 1, sizeof(seq_coor_t) );
    
    k_offset = (band_size >> 1) + 1;
    
    // We should probably use hashmap to store the backtracing information to save memory allocation time
    // This O(MN) block allocation scheme is convient for now but it is slower for very long sequences
    d_path = calloc( max_d * (band_size + 1 ) * 2 + 1, sizeof(d_path_data) );
    
    aln_path = calloc( q_len + t_len + 1, sizeof(path_point) );

    align_rtn = calloc( 1, sizeof(alignment));
    align_rtn->t_aln_str = calloc( q_len + t_len + 1, sizeof(char));
    align_rtn->q_aln_str = calloc( q_len + t_len + 1, sizeof(char));
    align_rtn->aln_str_size = 0;
    align_rtn->aln_q_s = 0;
    align_rtn->aln_q_e = 0;
    align_rtn->aln_t_s = 0;
    align_rtn->aln_t_e = 0;

    //printf("max_d: %lu, band_size: %lu\n", max_d, band_size);
    best_m = -1;
    min_k = 0;
    max_k = 0;
    
    for (d = 0; d < max_d; d ++ ) {
        max_k = max_k > max_d ? max_d : max_k;
        min_k = min_k < -max_d ? -max_d : min_k;

        if (max_k - min_k > 500) {
            //printf("band too wide\n");
            //fflush(stdout);
            break;
        }
        if (abs(max_k) > k_offset - 1  || abs(min_k) > k_offset - 1 ) break;

        //printf("d, min_k, max_k, x1, lim: %ld %ld %ld %ld %ld\n", d, min_k, max_k, band_size + max_k + k_offset, max_d * (band_size * 2 + 1) );
        //if (d * band_size + max_k + k_offset  >=  max_d * (band_size + 1) * 2 + 1 ) break;
        //if (d * band_size + min_k + k_offset  <  0 ) break;
 
        for (k = min_k; k <= max_k;  k += 2) {

            if ( k == min_k || k != max_k && V[ k - 1 + k_offset ] < V[ k + 1 + k_offset] ) {
                pre_k = k + 1;
                x = V[ k + 1 + k_offset];
            } else {
                pre_k = k - 1;
                x = V[ k - 1 + k_offset] + 1;
            }
            y = x - k;
            d_path[ d * band_size + k + k_offset ].x1 = x;
            d_path[ d * band_size + k + k_offset ].y1 = y;

            //printf("1 x,y,d,k = %ld %ld %ld %ld\n", x, y, d, k);
            //printf("%c %c\n", query_seq[x], target_seq[y]);

            while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){
                x++;
                y++;
            }

            //printf("2 x,y, q_len, t_len = %ld %ld %ld %ld\n\n", x, y, q_len, t_len);
            
            d_path[ d * band_size + k + k_offset ].x2 = x;
            d_path[ d * band_size + k + k_offset ].y2 = y;
            d_path[ d * band_size + k + k_offset ].pre_k = pre_k;

            V[ k + k_offset ] = x;
            U[ k + k_offset ] = x + y;
            
            if ( x + y > best_m) {
                best_m = x + y;
            }

            if ( x >= q_len || y >= t_len) {
                aligned = true;
                break;
            }
        }
        
        // For banding
        new_min_k = max_k;
        new_max_k = min_k;

        for (k2 = min_k; k2 <= max_k;  k2 += 2) {
            if (U[ k2 + k_offset] >= best_m - band_tolerance ) {
                if ( k2 < new_min_k ) {
                    new_min_k = k2;
                }
                if ( k2 > new_max_k ) {
                    new_max_k = k2;
                }
            }
        }
        
        max_k = new_max_k + 1;
        min_k = new_min_k - 1;
        
        // For no banding
        // max_k ++;
        // min_k --;

        // For debuging 
        //printf("min_max_k,d, %ld %ld %ld\n", min_k, max_k, d);
        
        if (aligned == true) {
            align_rtn->aln_q_e = x;
            align_rtn->aln_t_e = y;
            align_rtn->dist = d;
            align_rtn->aln_str_size = (x + y + d) / 2;
            align_rtn->aln_q_s = 0;
            align_rtn->aln_t_s = 0;

            if (get_aln_str > 0) {
                cd = d;
                ck = k;
                aln_path_idx = 0;
                while (cd >= 0 && aln_path_idx < q_len + t_len + 1) {    
                    aln_path[aln_path_idx].x = d_path[ cd * band_size + ck + k_offset ].x2;
                    aln_path[aln_path_idx].y = d_path[ cd * band_size + ck + k_offset ].y2;
                    aln_path_idx ++;
                    aln_path[aln_path_idx].x = d_path[ cd * band_size + ck + k_offset ].x1;
                    aln_path[aln_path_idx].y = d_path[ cd * band_size + ck + k_offset ].y1;
                    aln_path_idx ++;
                    ck = d_path[ cd * band_size + ck + k_offset ].pre_k;
                    cd -= 1;
                }
                aln_path_idx --;
                cx = aln_path[aln_path_idx].x;
                cy = aln_path[aln_path_idx].y;
                align_rtn->aln_q_s = cx;
                align_rtn->aln_t_s = cy;
                aln_pos = 0;
                while ( aln_path_idx > 0 ) {
                    aln_path_idx --;
                    nx = aln_path[aln_path_idx].x;
                    ny = aln_path[aln_path_idx].y;
                    if (cx == nx && cy == ny){
                        continue;
                    }
                    if (nx == cx && ny != cy){ //advance in y
                        for (i = 0; i <  ny - cy; i++) {
                            align_rtn->q_aln_str[aln_pos + i] = '-';
                        }
                        for (i = 0; i <  ny - cy; i++) {
                            align_rtn->t_aln_str[aln_pos + i] = target_seq[cy + i];
                        }
                        aln_pos += ny - cy;
                    } else if (nx != cx && ny == cy){ //advance in x
                        for (i = 0; i <  nx - cx; i++) {
                            align_rtn->q_aln_str[aln_pos + i] = query_seq[cx + i];
                        }
                        for (i = 0; i <  nx - cx; i++) {
                            align_rtn->t_aln_str[aln_pos + i] = '-';
                        }
                        aln_pos += nx - cx;
                    } else {
                        for (i = 0; i <  nx - cx; i++) {
                            align_rtn->q_aln_str[aln_pos + i] = query_seq[cx + i];
                        }
                        for (i = 0; i <  ny - cy; i++) {
                            align_rtn->t_aln_str[aln_pos + i] = target_seq[cy + i];
                        }
                        aln_pos += ny - cy;
                    }
                    cx = nx;
                    cy = ny;
                }
                align_rtn->aln_str_size = aln_pos;
            }
            break;
        }
    }

    free(V);
    free(U);
    free(d_path);
    free(aln_path);
    return align_rtn;
}


void free_alignment(alignment * aln) {
    free(aln->q_aln_str);
    free(aln->t_aln_str);
    free(aln);
}
