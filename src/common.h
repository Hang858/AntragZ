#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stddef.h>

#define ANTRAG_D 512        // n
#define ANTRAG_LOGD 9    // log2(n)
#define ANTRAG_Q 12289      // q
#define ANTRAG_ALPHA 1.17   // alpha
#define ANTRAG_ALPHAEPS 0.0
#define ANTRAG_XI 1.0 / 3.0      // xi, 通常接近1
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEPTH_INT_FG   4
#define fpr_add(a, b) ((a) + (b))
#define fpr_sub(a, b) ((a) - (b))
#define fpr_mul(a, b) ((a) * (b))
#define fpr_sqr(a)    ((a) * (a))
#define fpr_inv(a)    (1.0 / (a))
#define fpr_lt(a, b)  ((a) < (b))
#define fpr_rint(a)   rint(a)
#define fpr_of(a)     ((double)(a))
#define fpr_zero      0.0
#define fpr_one       1.0
#define fpr_two       2.0
#define fpr_onehalf   0.5
#define fpr_ptwo31    2147483648.0
#define fpr_ptwo31m1  2147483647.0
#define fpr_mtwo31m1  -2147483647.0
#define fpr_ptwo63m1  9223372036854775807.0
#define fpr_mtwo63m1  -9223372036854775807.0

#define MKN(logn) ((size_t)1 << (logn))

typedef struct {
    double coeffs[ANTRAG_D];
} poly;

typedef struct {
    int64_t *coeffs;
    size_t n; // 维度 (2^k)
} poly_int;

typedef struct {
    int8_t f[ANTRAG_D];  // 最终的整数 f
    int8_t g[ANTRAG_D];  // 最终的整数 g
    poly b10;            // 临时变量 / 频域 f
    poly b11;            // 临时变量 / 频域 g
} secret_key_fg;

typedef struct {
    uint32_t p;
    uint32_t g;
    uint32_t s;
} small_prime;



#endif