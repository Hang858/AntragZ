#ifndef MATH_ZINT_H
#define MATH_ZINT_H

#include <stdint.h>
#include <stddef.h>
#include "common.h"


uint32_t zint_sub(uint32_t *a, const uint32_t *b, size_t len, uint32_t ctl);
uint32_t zint_mul_small(uint32_t *m, size_t mlen, uint32_t x);
void zint_add_mul_small(uint32_t *x, const uint32_t *y, size_t len, uint32_t s);
void zint_norm_zero(uint32_t *x, const uint32_t *p, size_t len);
void zint_negate(uint32_t *a, size_t len, uint32_t ctl);


uint32_t zint_mod_small_unsigned(const uint32_t *d, size_t dlen, uint32_t p, uint32_t p0i, uint32_t R2);
uint32_t zint_mod_small_signed(const uint32_t *d, size_t dlen, uint32_t p, uint32_t p0i, uint32_t R2, uint32_t Rx);
void zint_rebuild_CRT(uint32_t *xx, size_t xlen, size_t xstride, size_t num, const small_prime *primes, int normalize_signed, uint32_t *tmp);

uint32_t zint_co_reduce(uint32_t *a, uint32_t *b, size_t len, int64_t xa, int64_t xb, int64_t ya, int64_t yb);
void zint_co_reduce_mod(uint32_t *a, uint32_t *b, const uint32_t *m, size_t len, uint32_t m0i, int64_t xa, int64_t xb, int64_t ya, int64_t yb);
void zint_finish_mod(uint32_t *a, size_t len, const uint32_t *m, uint32_t neg);
int zint_bezout(uint32_t *u, uint32_t *v, const uint32_t *x, const uint32_t *y, size_t len, uint32_t *tmp);


void zint_add_scaled_mul_small(uint32_t *x, size_t xlen, const uint32_t *y, size_t ylen, int32_t k, uint32_t sch, uint32_t scl);
void zint_sub_scaled(uint32_t *x, size_t xlen, const uint32_t *y, size_t ylen, uint32_t sch, uint32_t scl);

int32_t zint_one_to_plain(const uint32_t *x);
void poly_big_to_fp(double *d, const uint32_t *f, size_t flen, size_t fstride, unsigned logn);
int poly_big_to_small(int8_t *d, const uint32_t *s, int lim, unsigned logn);

#endif