#ifndef ENCODE_H
#define ENCODE_H
#include <stddef.h>
#include <stdint.h>

size_t compress_sig(void *out, size_t max_out_len, const int16_t *x);

size_t decompress_sig(int16_t *x, const uint8_t *in, size_t max_in_len);

void encode_public_key(uint8_t *out, const int16_t *h);

int decode_public_key(int16_t *h, const uint8_t *in, size_t len);

#endif