
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "common.h"
consensus_data * generate_consensus( char ** input_seq,
                           unsigned int n_seq,
                           unsigned min_cov,
                           unsigned K,
                           double min_idt);
void free_consensus_data( consensus_data * consensus );

void main() {
    unsigned int j;
    char small_buffer[1024];
    char big_buffer[65536];
    char ** input_seq;
    char ** seq_id;
    int seq_count;
    consensus_data * consensus;

    input_seq = calloc( 501, sizeof(char *));
    seq_id = calloc( 501, sizeof(char *));

    while(1) {
        seq_count = 0;
        while (1) {

            scanf("%s", small_buffer);
            seq_id[seq_count] = calloc( strlen(small_buffer) + 1, sizeof(char));
            strcpy(seq_id[seq_count], small_buffer);

            scanf("%s", big_buffer);
            input_seq[seq_count] = calloc( strlen(big_buffer) + 1 , sizeof(char));
            strcpy(input_seq[seq_count], big_buffer);

            if (strcmp(seq_id[seq_count], "+") == 0) {
                break;
            }
            if (strcmp(seq_id[seq_count], "-") == 0) {
                break;
            }
            //printf("%s\n", seq_id[seq_count]);
            seq_count += 1;
            if (seq_count > 500) break;
        }
        //printf("sc: %d\n", seq_count);
        if (seq_count < 10 && strcmp(seq_id[seq_count], "-") != 0 ) continue;
        if (seq_count < 10 && strcmp(seq_id[seq_count], "-") == 0 ) break;

        consensus = generate_consensus(input_seq, seq_count, 8, 8, 0.8);

        if (strlen(consensus->sequence) > 500) {
            printf(">%s\n%s\n", seq_id[0], consensus->sequence);
        }
        fflush(stdout);
	free_consensus_data(consensus);
        for (j=0; j < seq_count; j++) {
            free(seq_id[j]);
            free(input_seq[j]);
        };

    }
    for (j=0; j < seq_count; j++) {
        free(seq_id[j]);
        free(input_seq[j]);
    };
    free(seq_id);
    free(input_seq);
}
