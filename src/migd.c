#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "migd.h"

static int128_t calc_g_norm_sq_migd(int64_t b, int k) {
    int128_t sum = 0;
    int128_t bj = 1;
    for (int i = 0; i < k; i++) {
        sum += bj * bj;
        if (i < k - 1) bj *= b;
    }
    return sum;
}


static void poly_adj(const int64_t *src, int64_t *dst, int n) {
    dst[0] = src[0];
    for (int i = 1; i < n; i++) {
        dst[i] = -src[n - i];
    }
}


static void poly_mul_add_accum(int64_t *target, const int64_t *a, const int64_t *b, int n) {

    int128_t *tmp = (int128_t *)calloc(n, sizeof(int128_t));
    if (!tmp) return;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int k = i + j;
            int128_t val = (int128_t)a[i] * (int128_t)b[j];
            
            if (k < n) tmp[k] += val;
            else       tmp[k - n] -= val;
        }
    }
    
    for(int i=0; i<n; i++) {
        target[i] += (int64_t)tmp[i];
    }
    free(tmp);
}

static void poly_decompose_gadget(const int64_t *poly_in, int64_t b, int k, int n, int64_t *decomposed_polys) {
    memset(decomposed_polys, 0, k * n * sizeof(int64_t));
    
    for (int i = 0; i < n; i++) {
        int64_t val = poly_in[i];
        
        for (int j = 0; j < k; j++) {
            int64_t r = val % b;
            if (r < 0) r += b;
            if (r >= b / 2) r -= b;
            
            decomposed_polys[j * n + i] = r;
            
            val = (val - r) / b;
        }
    }
}


int Run_MIGD(
    const int64_t *sigma_11, 
    const int64_t *sigma_12, 
    const int64_t *sigma_22,
    int128_t d,
    int64_t b, 
    MIGD_Output *out_A
) {

    memset(&out_A->x_ctx, 0, sizeof(MemContext));
    memset(&out_A->y_ctx, 0, sizeof(MemContext));
    
    int n = N_MAX;
    int k = K_VAL;

    printf("[MIGD] Starting Decomposition...\n");


    int128_t g_norm_sq = calc_g_norm_sq_migd(b, k);
    int128_t target_d = d - g_norm_sq;
    

    if (target_d <= 0) {
        printf("[MIGD] Error: d is too small (d' <= 0).\n");
        return 0;
    }

 
    printf("[MIGD] Decomposing Sigma_11 (x)...\n");
    if (!EIGD_recursive_opt(sigma_11, n, L_VAL, 1, target_d, b, k, &out_A->x_ctx, 0)) {
        printf("[MIGD] EIGD for x failed.\n");
        return 0;
    }

    printf("[MIGD] Decomposing Sigma_12 (c)...\n");
    int64_t *neg_sigma_12 = (int64_t *)malloc(n * sizeof(int64_t));
    if (!neg_sigma_12) return 0;

    for(int i=0; i<n; i++) neg_sigma_12[i] = -sigma_12[i];

    poly_decompose_gadget(neg_sigma_12, b, k, n, out_A->c_coeffs);
    free(neg_sigma_12);

    printf("[MIGD] Computing Schur Complement Pi...\n");
    
    int64_t *Pi = (int64_t *)calloc(n, sizeof(int64_t));
    int64_t *tmp_adj = (int64_t *)malloc(n * sizeof(int64_t));
    
    if (!Pi || !tmp_adj) {
        free(Pi); free(tmp_adj);
        return 0;
    }

    for(int i=0; i<n; i++) Pi[i] = sigma_22[i];


    for (int j = 0; j < k; j++) {

        int64_t *c_j = out_A->c_coeffs + (j * n);
        poly_adj(c_j, tmp_adj, n);
        poly_mul_add_accum(Pi, tmp_adj, c_j, n);
    }
    
    free(tmp_adj);


    printf("[MIGD] Decomposing Pi (y)...\n");
    if (!EIGD_recursive_opt(Pi, n, L_VAL, 1, target_d, b, k, &out_A->y_ctx, 0)) {
        printf("[MIGD] EIGD for y failed (Schur complement invalid?).\n");
        free(Pi);
        return 0;
    }

    free(Pi);
    printf("[MIGD] Success! Matrix A constructed.\n");
    
    return 1;
}