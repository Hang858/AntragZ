#ifndef EIGD_H
#define EIGD_H

#include <stdint.h>
#include "common.h"



int EIGD_recursive_opt(const int64_t *input_f, int n, int current_l, int stride, 
                       int128_t d, int64_t b, int k, MemContext *ctx, int stack_offset);

void init_compression_tables(void);

#endif // EIGD_H