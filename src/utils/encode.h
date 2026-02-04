#ifndef ENCODE_H
#define ENCODE_H
#include <stddef.h>
#include <stdint.h>

size_t compress_sig(void *out, size_t max_out_len, const int16_t *x);

#endif