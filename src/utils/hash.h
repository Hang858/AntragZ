#ifndef HASH_H
#define HASH_H
#include <stdint.h>
#include <stddef.h>
typedef struct {
    uint64_t ctx[26];
} shake256incctx;

void hash_to_point(int16_t *x, const uint8_t *r, size_t r_len, const uint8_t *msg, size_t msg_len);

#endif