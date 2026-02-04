#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "prematrix.h"
#include "fft.h"
#include "common.h"
#include "migd.h"
#include "poly.h"
#define PREMATRIX_N 512


void run_cholesky(poly *a11, poly *a21, poly *a22,
                  const poly *p11, const poly *p21, const poly *p22);

int Run_PreMatrix(
    const int8_t *f, const int8_t *g, 
    const int8_t *F, const int8_t *G, 
    PreMatrix_Output *out
) {
    int n = ANTRAG_D;
    
    // 安全检查：头文件不匹配会导致内存破坏
    if (n > PREMATRIX_N) {
        printf("[ERROR] ANTRAG_D > PREMATRIX_N. Recompile with updated prematrix.h\n");
        return 0;
    }

    int64_t p = 1LL << 28; 
    
    // FFT
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

    poly *num = malloc(sizeof(poly));
    poly_add_muladj_fft(num, tp3, tp4, tp1, tp2, ANTRAG_LOGD); 

    poly *den = malloc(sizeof(poly));
    poly_mulselfadj_fft(tp1, ANTRAG_LOGD); 
    poly_mulselfadj_fft(tp2, ANTRAG_LOGD); 
    
    for(int i=0; i<n; i++) den->coeffs[i] = tp1->coeffs[i] + tp2->coeffs[i];

    for(int i=0; i < n/2; i++) {
        double nr = num->coeffs[i];
        double ni = num->coeffs[i + n/2];
        double dr = den->coeffs[i]; 
        tp1->coeffs[i] = nr / dr;
        tp1->coeffs[i + n/2] = ni / dr;
    }
    
    invFFT(tp1, ANTRAG_LOGD);
    
    for(int i=0; i<n; i++) {
        out->u_hat_num[i] = (int64_t)round(tp1->coeffs[i] * p);
    }
    
    free(num); free(den); 
    free(tp1); free(tp2); free(tp3); free(tp4);

    // 计算 Sigma
    int64_t *R1 = malloc(n * sizeof(int64_t));
    int64_t *R2 = malloc(n * sizeof(int64_t));
    int64_t *f_int = malloc(n * sizeof(int64_t));
    int64_t *g_int = malloc(n * sizeof(int64_t));
    int128_t *tmp128 = calloc(n, sizeof(int128_t));
    
    for(int i=0; i<n; i++) { f_int[i] = f[i]; g_int[i] = g[i]; }

    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_acc_128(tmp128, out->u_hat_num, f_int, n);
    for(int i=0; i<n; i++) R1[i] = (int64_t)F[i] * p - (int64_t)tmp128[i];
    
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_acc_128(tmp128, out->u_hat_num, g_int, n);
    for(int i=0; i<n; i++) R2[i] = (int64_t)G[i] * p - (int64_t)tmp128[i];

    int128_t *S11 = calloc(n, sizeof(int128_t));
    int128_t *S12 = calloc(n, sizeof(int128_t));
    int128_t *S22 = calloc(n, sizeof(int128_t));
    int128_t p2 = (int128_t)p * p;
    
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_adj_acc_128(tmp128, f_int, f_int, n);
    for(int i=0; i<n; i++) S11[i] += tmp128[i] * p2;
    poly_mul_adj_acc_128(S11, R1, R1, n);
    
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_adj_acc_128(tmp128, g_int, g_int, n);
    for(int i=0; i<n; i++) S22[i] += tmp128[i] * p2;
    poly_mul_adj_acc_128(S22, R2, R2, n);
    
    memset(tmp128, 0, n * sizeof(int128_t));
    poly_mul_adj_acc_128(tmp128, f_int, g_int, n);
    for(int i=0; i<n; i++) S12[i] += tmp128[i] * p2;
    poly_mul_adj_acc_128(S12, R1, R2, n);

    free(R1); free(R2); free(f_int); free(g_int);

    // 计算 P
    int64_t s0 = 131;
    int64_t B_val = (int64_t)p * (s0 - 1);
    int128_t B_sq = (int128_t)B_val * B_val; // 注意这里是 (p*(s0-1))^2
    
    for(int i=0; i<n; i++) {
        S11[i] = -S11[i];
        if (i == 0) S11[0] += B_sq;
    }
    for(int i=0; i<n; i++) {
        S22[i] = -S22[i];
        if (i == 0) S22[0] += B_sq;
    }
    
    // P21 = -Sigma12* = -S12*
    // Code check: S12* [k] is S12[n-k] (with sign flip logic). 
    // Wait. S12[k] is coefficient of x^k.
    // (S12)^*[k] = -S12[n-k]. (for k!=0)
    // So P21[k] = -(-S12[n-k]) = S12[n-k].
    // This is mathematically correct.
    int128_t *P21 = malloc(n * sizeof(int128_t));
    P21[0] = -S12[0];
    for(int i=1; i<n; i++) {
        P21[i] = S12[n - i]; 
    }
    free(S12); 

    // Cholesky
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
        printf("[ERROR] Cholesky produced NaN. Matrix not PD.\n");
        // ... (释放内存)
        free(poly_p11); free(poly_p21); free(poly_p22);
        free(poly_c11); free(poly_c21); free(poly_c22);
        free(S11); free(S22); free(P21); free(tmp128);
        return 0;
    }
    
    // Check for NaN
    if (isnan(poly_c11->coeffs[0])) {
        printf("[ERROR] Cholesky produced NaN. Matrix not positive definite.\n");
        // Clean up and fail
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

    // Delta
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
    
    // MIGD Inputs
    int64_t *d11_64 = malloc(n * sizeof(int64_t));
    int64_t *d12_64 = malloc(n * sizeof(int64_t)); 
    int64_t *d22_64 = malloc(n * sizeof(int64_t));
    
    // d12 = Delta_12 = Delta_21^*
    // P21 holds Delta_21
    // Delta_21^*[k] is P21[0] or -P21[n-k]
    for(int i=0; i<n; i++) {
        d11_64[i] = (int64_t)S11[i];
        d22_64[i] = (int64_t)S22[i];
        if (i==0) d12_64[0] = (int64_t)P21[0];
        else      d12_64[i] = -(int64_t)P21[n-i];
    }
    
    free(S11); free(S22); free(P21); free(tmp128);

    int128_t d_migd = (int128_t)2 * p * (p * (s0 - 1)) - 1;
    int64_t b_param = 32768; 

    int res = Run_MIGD(d11_64, d12_64, d22_64, d_migd, b_param, &out->migd_key);
    out->valid = res;

    free(d11_64); free(d12_64); free(d22_64);
    
    return res;
}