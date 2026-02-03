#ifndef MATH_POLY_H
#define MATH_POLY_H

#include <stdint.h>
#include "common.h" // 获取 ANTRAG_D, ANTRAG_Q 定义

// 模 q 归约
static inline int32_t mod_q(int32_t a) {
    int32_t r = a % ANTRAG_Q;
    return r < 0 ? r + ANTRAG_Q : r;
}

// 整数求逆 (扩展欧几里得)
int16_t mod_inverse(int16_t a, int16_t m);

// 多项式乘法: res = a * b mod (x^N + 1, q)
void poly_mul_mod(int16_t *res, const int16_t *a, const int16_t *b);

// 多项式求逆: out = f^-1 mod (x^N + 1, q)
// 返回 1 成功，0 失败
int poly_inv_mod_q(int16_t *out, const int8_t *f_in);

#endif // MATH_POLY_H