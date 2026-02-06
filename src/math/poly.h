#ifndef MATH_POLY_H
#define MATH_POLY_H

#include <stdint.h>
#include "common.h"

static inline int32_t mod_q(int32_t a) {
    int32_t r = a % ANTRAG_Q;
    return r < 0 ? r + ANTRAG_Q : r;
}

int16_t mod_inverse(int16_t a, int16_t m);

void poly_mul_mod(int16_t *res, const int16_t *a, const int16_t *b);

int poly_inv_mod_q(int16_t *out, const int8_t *f_in);


void poly_adj(const int64_t *src, int64_t *dst, int n);

void poly_decompose_gadget(const int64_t *poly_in, int64_t b, int k, int n, int64_t *decomposed_polys);


void poly_mul_acc_64(int64_t *target, const int64_t *a, const int64_t *b, int n);

void poly_mul_acc_128(int128_t *res, const int64_t *a, const int64_t *b, int n);


void poly_mul_adj_acc_128(int128_t *res, const int64_t *a, const int64_t *b, int n);

void poly_mul_sub_scaled_128(int128_t *res, const int64_t *a, const int64_t *b, int64_t scale, int n);
void poly_mul_int8_int64_to_128_acc(int128_t *res, const int8_t *a, const int64_t *b, int n);
void poly_mul_int8_int64_acc(int64_t *res, const int8_t *a, const int64_t *b, int n);
void int128_to_poly_double(poly *out, const int128_t *in, int n);
void poly_double_to_int64(int64_t *out, const poly *in, int n);

#endif