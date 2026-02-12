#ifndef POLY_NTT_512_H
#define POLY_NTT_512_H

#include <stdint.h>

// Parameters
#define POLY_512_N 512
#define POLY_512_Q 12289

/**
 * @brief Initialize hardware context for 512-NTT.
 * Must be called once before using poly_ntt_512_mul.
 */
void poly_ntt_512_init(void);

/**
 * @brief Perform 512-dimension polynomial multiplication.
 * res = a * b mod (X^512 + 1) mod 12289
 * * @param res Output array (512 int16_t)
 * @param a   Input array A (512 int16_t)
 * @param b   Input array B (512 int16_t)
 * @return    0 on success, -1 on failure
 */
int poly_ntt_512_mul(int16_t *res, const int16_t *a, const int16_t *b);

#endif // POLY_NTT_512_H