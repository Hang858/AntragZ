#ifndef EIGD_OPT_H
#define EIGD_OPT_H

#include <stdint.h>

#define N_MAX 512
#define K_VAL 3
#define L_VAL 9
#define OUTPUT_ROWS 28 // 3*(9-1) + 4


#define COMPRESSED_POOL_SIZE 1024 


#define STACK_BUF_SIZE 512
#define WORK_BUF_SIZE 1600


typedef __int128_t int128_t;

typedef struct {

    int16_t coeffs[COMPRESSED_POOL_SIZE];
    
    int64_t base_case[4];

    int64_t input_stack[STACK_BUF_SIZE];
    int64_t workspace[WORK_BUF_SIZE];
} MemContext;


int EIGD_recursive_opt(const int64_t *input_f, int n, int current_l, int stride, 
                       int128_t d, int64_t b, int k, MemContext *ctx, int stack_offset);

void init_compression_tables(void);

#endif // EIGD_OPT_H