#ifndef EIGD_H
#define EIGD_H

#include <stdint.h>
#include "common.h"



int EIGD_recursive_opt(const int64_t *input_f, int n, int current_l, int stride, 
                       int128_t d, int64_t b, int k, MemContext *ctx, int stack_offset);

int64_t matrix_get(MemContext *ctx, int row, int col);

#endif // EIGD_H