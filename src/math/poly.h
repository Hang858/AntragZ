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


/**
 * @brief 计算多项式的共轭 f*(x) = f(1/x)
 * dst[0] = src[0], dst[i] = -src[n-i]
 */
void poly_adj(const int64_t *src, int64_t *dst, int n);

/**
 * @brief Gadget 分解 (用于 MIGD)
 * 将 poly_in 分解为 k 个多项式，每个系数在 [-b/2, b/2)
 */
void poly_decompose_gadget(const int64_t *poly_in, int64_t b, int k, int n, int64_t *decomposed_polys);

// ==========================================
// 乘法与累加 (Multiplication & Accumulation)
// ==========================================

/**
 * @brief 乘法累加: target += a * b
 * @note 内部使用 int128 防止溢出，结果存回 int64
 * (原 migd.c 版本，较安全)
 */
void poly_mul_acc_64(int64_t *target, const int64_t *a, const int64_t *b, int n);

/**
 * @brief 高精度乘法累加: res += a * b
 * @note 结果存入 int128 数组 (原 prematrix.c 版本)
 */
void poly_mul_acc_128(int128_t *res, const int64_t *a, const int64_t *b, int n);

/**
 * @brief 高精度共轭乘法累加: res += a * adj(b)
 * @note 结果存入 int128 数组
 */
void poly_mul_adj_acc_128(int128_t *res, const int64_t *a, const int64_t *b, int n);

void poly_mul_sub_scaled_128(int128_t *res, const int64_t *a, const int64_t *b, int64_t scale, int n);
void poly_mul_int8_int64_to_128_acc(int128_t *res, const int8_t *a, const int64_t *b, int n);
void poly_mul_int8_int64_acc(int64_t *res, const int8_t *a, const int64_t *b, int n);

// ==========================================
// 类型转换 (Conversion)
// ==========================================

void int128_to_poly_double(poly *out, const int128_t *in, int n);
void poly_double_to_int64(int64_t *out, const poly *in, int n);

#endif // MATH_POLY_H