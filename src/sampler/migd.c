#include <stdlib.h>
#include <string.h>
#include "migd.h"
#include "eigd.h"
#include "poly.h"

#ifdef ZITAKA_DEBUG
static int64_t get_max_abs_64(const int64_t *v, int n) {
    int64_t max = 0;
    for(int i=0; i<n; i++) {
        int64_t val = v[i] < 0 ? -v[i] : v[i];
        if(val > max) max = val;
    }
    return max;
}
#endif


static int128_t calc_g_norm_sq_migd(int64_t b, int k) {
    int128_t sum = 0;
    int128_t bj = 1;
    for (int i = 0; i < k; i++) {
        sum += bj * bj;
        if (i < k - 1) bj *= b;
    }
    return sum;
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


    int128_t g_norm_sq = calc_g_norm_sq_migd(b, k);
    int128_t target_d = d - g_norm_sq;

    ZITAKA_LOG("[MIGD] Start. d=%.1e, ||g||^2=%.1e, target_d=%.1e", 
               (double)d, (double)g_norm_sq, (double)target_d);
    
    
    // 对 sigma_11 进行 EIGD 分解
    // x * x* = target_d - Sigma_11
#ifdef ZITAKA_DEBUG
    int64_t max_s11 = get_max_abs_64(sigma_11, n);
    ZITAKA_LOG("[MIGD] EIGD(x) input max coeff: %ld", max_s11);
#endif

    if (!EIGD_recursive_opt(sigma_11, n, L_VAL, 1, target_d, b, k, &out_A->x_ctx, 0)) {
        ZITAKA_LOG("[MIGD] EIGD for x failed.\n");
        return 0;
    }

    // Gadget 分解求解 c_j
    int64_t *neg_sigma_12 = (int64_t *)malloc(n * sizeof(int64_t));
    if (!neg_sigma_12) return 0;

    for(int i=0; i<n; i++) neg_sigma_12[i] = -sigma_12[i];

    // sum (b^j * c_j) = -Sigma_12
    poly_decompose_gadget(neg_sigma_12, b, k, n, out_A->c_coeffs);
    free(neg_sigma_12);

    // 计算 Pi = Sigma_22 + sum(c_j* * c_j)
    int64_t *Pi = (int64_t *)calloc(n, sizeof(int64_t));
    int64_t *tmp_adj = (int64_t *)malloc(n * sizeof(int64_t));
    
    if (!Pi || !tmp_adj) {
        free(Pi); free(tmp_adj);
        return 0;
    }
    // 初始化 Pi = Sigma_22
    for(int i=0; i<n; i++) Pi[i] = sigma_22[i];

    // 累加 sum(c_j* * c_j)
    for (int j = 0; j < k; j++) {

        int64_t *c_j = out_A->c_coeffs + (j * n);
        poly_adj(c_j, tmp_adj, n);
        poly_mul_acc_64(Pi, tmp_adj, c_j, n);
    }
    
    free(tmp_adj);

#ifdef ZITAKA_DEBUG
    int64_t max_pi = get_max_abs_64(Pi, n);
    int128_t remaining_budget = target_d - (int128_t)Pi[0];
    ZITAKA_LOG("[MIGD] Pi computed. Max coeff: %ld. Pi[0]: %ld", max_pi, Pi[0]);
    
    if (remaining_budget < 0) {
        ZITAKA_LOG("[MIGD] Critical Warning: Pi[0] > target_d! (Diff: %.1e). EIGD(y) will fail.", (double)remaining_budget);
    }

#endif

    if (!EIGD_recursive_opt(Pi, n, L_VAL, 1, target_d, b, k, &out_A->y_ctx, 0)) {
        ZITAKA_LOG("[MIGD] Error: EIGD for y failed.");
        free(Pi);
        return 0;
    }

    free(Pi);
    ZITAKA_LOG("[MIGD] Success.");
    return 1;
}