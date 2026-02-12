#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "prematrix.h"
#include "fft.h"
#include "common.h"
#include "migd.h"
#include "poly.h"
#define PREMATRIX_N 512

#ifdef ZITAKA_DEBUG
static int64_t get_max_abs_128(const int128_t *v, int n) {
    int64_t max = 0;
    for(int i=0; i<n; i++) {
        int128_t val = v[i] < 0 ? -v[i] : v[i];
        if(val > (int128_t)INT64_MAX) val = INT64_MAX; // 截断显示
        if((int64_t)val > max) max = (int64_t)val;
    }
    return max;
}

#endif

void run_cholesky(poly *a11, poly *a21, poly *a22,
                  const poly *p11, const poly *p21, const poly *p22);

int Run_PreMatrix(
    const int8_t *f, const int8_t *g, 
    const int8_t *F, const int8_t *G, 
    PreMatrix_Output *out
) {
    int n = ANTRAG_D;
    ZITAKA_LOG("[PreMatrix] Start. n=%d", n);
    if (n > PREMATRIX_N) {
        ZITAKA_LOG("[ERROR] ANTRAG_D > PREMATRIX_N. Recompile required.");
        return 0;
    }

    int64_t p = 1LL << 28; 
    
    // FFT 计算 u_hat
    poly *tp1 = malloc(sizeof(poly));
    poly *tp2 = malloc(sizeof(poly));
    poly *tp3 = malloc(sizeof(poly));
    poly *tp4 = malloc(sizeof(poly));
    
    if(!tp1 || !tp2 || !tp3 || !tp4) return 0; 

    for(int i=0; i<n; i++) {
        tp1->coeffs[i] = (double)f[i];
        tp2->coeffs[i] = (double)g[i];
        tp3->coeffs[i] = (double)F[i];
        tp4->coeffs[i] = (double)G[i];
    }

    FFT(tp1, ANTRAG_LOGD); FFT(tp2, ANTRAG_LOGD);
    FFT(tp3, ANTRAG_LOGD); FFT(tp4, ANTRAG_LOGD);
    // 计算 u 的分子
    poly *num = malloc(sizeof(poly));
    poly_add_muladj_fft(num, tp3, tp4, tp1, tp2, ANTRAG_LOGD); 
    // 计算 u 的分母
    poly *den = malloc(sizeof(poly));
    poly_mulselfadj_fft(tp1, ANTRAG_LOGD); 
    poly_mulselfadj_fft(tp2, ANTRAG_LOGD); 
    
    for(int i=0; i<n; i++) den->coeffs[i] = tp1->coeffs[i] + tp2->coeffs[i];
    // 复数除法，算u
    for(int i=0; i < n/2; i++) {
        double nr = num->coeffs[i];
        double ni = num->coeffs[i + n/2];
        double dr = den->coeffs[i]; 
        tp1->coeffs[i] = nr / dr;
        tp1->coeffs[i + n/2] = ni / dr;
    }
    
    invFFT(tp1, ANTRAG_LOGD);

    // 存储 u_hat_num
    for(int i=0; i<n; i++) {
        out->u_hat_num[i] = (int64_t)round(tp1->coeffs[i] * p);
    }

#ifdef ZITAKA_DEBUG
    int64_t max_u = 0;
    for(int i=0; i<n; i++) {
        int64_t val = out->u_hat_num[i] < 0 ? -out->u_hat_num[i] : out->u_hat_num[i];
        if(val > max_u) max_u = val;
    }
    ZITAKA_LOG("[PreMatrix] u_hat computed. Max coeff: %ld", max_u);
#endif
    
    free(num); free(den); 
    free(tp1); free(tp2); free(tp3); free(tp4);

    // 计算 Sigma (S11, S12, S22)
    int64_t *R1 = malloc(n * sizeof(int64_t));
    int64_t *R2 = malloc(n * sizeof(int64_t));
    int64_t *f_int = malloc(n * sizeof(int64_t));
    int64_t *g_int = malloc(n * sizeof(int64_t));
    int128_t *tmp128 = calloc(n, sizeof(int128_t));
    
    for(int i=0; i<n; i++) { f_int[i] = f[i]; g_int[i] = g[i]; }

    //计算 p*F - u_hat_num * f
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_acc_128(tmp128, out->u_hat_num, f_int, n);
    for(int i=0; i<n; i++) R1[i] = (int64_t)F[i] * p - (int64_t)tmp128[i];
    //计算 G*p - u_hat * g
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_acc_128(tmp128, out->u_hat_num, g_int, n);
    for(int i=0; i<n; i++) R2[i] = (int64_t)G[i] * p - (int64_t)tmp128[i];

    int128_t *S11 = calloc(n, sizeof(int128_t));
    int128_t *S12 = calloc(n, sizeof(int128_t));
    int128_t *S22 = calloc(n, sizeof(int128_t));
    int128_t p2 = (int128_t)p * p;
    // 计算 Sigma = p^2 * B * B^T
    // 计算 p^2 * f * f^*
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_adj_acc_128(tmp128, f_int, f_int, n);
    for(int i=0; i<n; i++) S11[i] += tmp128[i] * p2;
    // 累加 R1 * R1^*
    poly_mul_adj_acc_128(S11, R1, R1, n);
    // 计算 p^2 * g * g^*
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_adj_acc_128(tmp128, g_int, g_int, n);
    for(int i=0; i<n; i++) S22[i] += tmp128[i] * p2;
    // 累加 R2 * R2^*
    poly_mul_adj_acc_128(S22, R2, R2, n);
    // 同上，计算矩阵 sigma 的各个元素
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_adj_acc_128(tmp128, f_int, g_int, n);
    for(int i=0; i<n; i++) S12[i] += tmp128[i] * p2;
    poly_mul_adj_acc_128(S12, R1, R2, n);

    free(R1); free(R2); free(f_int); free(g_int);

    // 构造 Cholesky 输入矩阵 P = B^2*I - Sigma
    int64_t s0 = 131;
    int64_t B_val = (int64_t)p * (s0 - 1);
    int128_t B_sq = (int128_t)B_val * B_val;
    
    for(int i=0; i<n; i++) {
        S11[i] = -S11[i];
        if (i == 0) S11[0] += B_sq;
    }
    for(int i=0; i<n; i++) {
        S22[i] = -S22[i];
        if (i == 0) S22[0] += B_sq;
    }
    
    int128_t *P21 = malloc(n * sizeof(int128_t));
    P21[0] = -S12[0];
    for(int i=1; i<n; i++) {
        P21[i] = S12[n - i]; 
    }
    free(S12);

#ifdef ZITAKA_DEBUG
    ZITAKA_LOG("[PreMatrix] Cholesky input P prepared.");
    ZITAKA_LOG("  Max(|P11|): %ld", get_max_abs_128(S11, n));
    ZITAKA_LOG("  Max(|P21|): %ld", get_max_abs_128(P21, n));
    ZITAKA_LOG("  Max(|P22|): %ld", get_max_abs_128(S22, n));
#endif

    // 执行 Cholesky 分解
    poly *poly_p11 = malloc(sizeof(poly));
    poly *poly_p21 = malloc(sizeof(poly));
    poly *poly_p22 = malloc(sizeof(poly));
    
    int128_to_poly_double(poly_p11, S11, n);
    int128_to_poly_double(poly_p21, P21, n);
    int128_to_poly_double(poly_p22, S22, n);
    
    poly *poly_c11 = malloc(sizeof(poly));
    poly *poly_c21 = malloc(sizeof(poly));
    poly *poly_c22 = malloc(sizeof(poly));

    run_cholesky(poly_c11, poly_c21, poly_c22, poly_p11, poly_p21, poly_p22);

    if (isnan(poly_c11->coeffs[0]) || isnan(poly_c22->coeffs[0])) {
        ZITAKA_LOG("[ERROR] Cholesky produced NaN. Matrix not Positive Definite.");

        free(poly_p11); free(poly_p21); free(poly_p22);
        free(poly_c11); free(poly_c21); free(poly_c22);
        free(S11); free(S22); free(P21); free(tmp128);
        return 0;
    }

    poly_double_to_int64(out->c11, poly_c11, n);
    poly_double_to_int64(out->c21, poly_c21, n);
    poly_double_to_int64(out->c22, poly_c22, n);

    free(poly_p11); free(poly_p21); free(poly_p22);
    free(poly_c11); free(poly_c21); free(poly_c22);

    // 计算 Delta = C*C^T - P
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_adj_acc_128(tmp128, out->c11, out->c11, n);
    for(int i=0; i<n; i++) S11[i] = tmp128[i] - S11[i]; 
    
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_adj_acc_128(tmp128, out->c21, out->c21, n);
    poly_mul_adj_acc_128(tmp128, out->c22, out->c22, n);
    for(int i=0; i<n; i++) S22[i] = tmp128[i] - S22[i]; 
    
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_adj_acc_128(tmp128, out->c21, out->c11, n);
    for(int i=0; i<n; i++) P21[i] = tmp128[i] - P21[i]; 
    
    int64_t *d11_64 = malloc(n * sizeof(int64_t));
    int64_t *d12_64 = malloc(n * sizeof(int64_t)); 
    int64_t *d22_64 = malloc(n * sizeof(int64_t));
    
    for(int i=0; i<n; i++) {
        d11_64[i] = (int64_t)S11[i];
        d22_64[i] = (int64_t)S22[i];
        if (i==0) d12_64[0] = (int64_t)P21[0];
        else      d12_64[i] = -(int64_t)P21[n-i];
    }
#ifdef ZITAKA_DEBUG
    ZITAKA_LOG("[PreMatrix] Delta (MIGD input) ready.");
    ZITAKA_LOG("  Max(|D11|): %ld", get_max_abs_128(S11, n));
    ZITAKA_LOG("  Max(|D12|): %ld", get_max_abs_128(P21, n));
    ZITAKA_LOG("  Max(|D22|): %ld", get_max_abs_128(S22, n));
#endif
    
    free(S11); free(S22); free(P21); free(tmp128);

    // 计算 d_migd
    int128_t d_migd = (int128_t)2 * p * (p * (s0 - 1)) - 1;
    int64_t b_param = 32768;
    ZITAKA_LOG("[PreMatrix] Calling Run_MIGD with d_migd=%.1e", (double)d_migd);

    int res = Run_MIGD(d11_64, d12_64, d22_64, d_migd, b_param, &out->migd_key);
    out->valid = res;

    free(d11_64); free(d12_64); free(d22_64);

    if(res) {
        ZITAKA_LOG("[PreMatrix] SUCCESS.");
    } else {
        ZITAKA_LOG("[PreMatrix] FAILED at MIGD stage.");
    }
    
    return res;
}