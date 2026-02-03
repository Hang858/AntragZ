#ifndef MIGD_H
#define MIGD_H

#include "eigd.h"

typedef struct {
    MemContext x_ctx;
    MemContext y_ctx;
    int64_t c_coeffs[K_VAL * N_MAX]; 
} MIGD_Output;

/**
 * @brief
 * * @param sigma_11
 * @param sigma_12
 * @param sigma_22
 * @param d
 * @param b
 * @param out_A
 * @return int
 */
int Run_MIGD(
    const int64_t *sigma_11, 
    const int64_t *sigma_12, 
    const int64_t *sigma_22,
    int128_t d,
    int64_t b, 
    MIGD_Output *out_A
);

#endif // MIGD_H