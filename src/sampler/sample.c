#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "common.h"
#include "zalcon_samp.h"
#include "migd.h"
#include "prematrix.h"
#include "poly.h"
#include "rng.h"
#include "eigd.h"
#include <stdio.h>
#include <math.h>

#ifdef ZITAKA_DEBUG

static void int128_to_str(int128_t n, char *buf) {
    if (n == 0) { strcpy(buf, "0"); return; }
    
    char temp[50];
    int i = 0;
    int sign = 0;
    if (n < 0) { sign = 1; n = -n; }
    
    // 提取每一位
    while (n > 0) {
        temp[i++] = (char)((n % 10) + '0');
        n /= 10;
    }
    
    int j = 0;
    if (sign) buf[j++] = '-';
    while (i > 0) {
        buf[j++] = temp[--i];
    }
    buf[j] = '\0';
}

static int128_t get_real_max_abs_128(const int128_t *v, int n) {
    int128_t max = 0;
    for(int i=0; i<n; i++) {
        int128_t val = v[i] < 0 ? -v[i] : v[i];
        if(val > max) max = val;
    }
    return max;
}

#endif

static void get_eigd_poly_from_ctx(const MemContext *ctx, int row, int64_t *out_poly, int n) {
    for(int i=0; i<n; i++) {
        out_poly[i] = matrix_get(ctx, row, i);
    }
}

void OfflineSamp(const PreMatrix_Output *key, int64_t *out_p1, int64_t *out_p2) {
    int n = ANTRAG_D;

    // v1, v2 是 128 位累加器，用于存储矩阵乘法结果 A * x
    int128_t *v1 = calloc(n, sizeof(int128_t));
    int128_t *v2 = calloc(n, sizeof(int128_t));
    
    int64_t *temp_x = malloc(n * sizeof(int64_t));     // 临时的高斯噪声向量 x
    int64_t *temp_poly = malloc(n * sizeof(int64_t));
    
    if (!v1 || !v2 || !temp_x || !temp_poly) goto cleanup;

    for(int i=0; i<n; i++) v1[i] += SampleLW(); // 采样 标准差为Lr0
    
    for(int i=0; i<n; i++) v2[i] += SampleLW(); 

    for(int i=0; i<n; i++) temp_x[i] = SampleLW();
    // v1 += C11 * x1
    poly_mul_acc_128(v1, key->c11, temp_x, n);
    // v2 += C21 * x1
    poly_mul_acc_128(v2, key->c21, temp_x, n);
    

    for(int i=0; i<n; i++) temp_x[i] = SampleLW();
    // v2 += C22 * x2
    poly_mul_acc_128(v2, key->c22, temp_x, n);


    int64_t b_pow = 1;
    const int64_t b = 32768;
    
    for(int j=0; j<K_VAL; j++) {

        for(int i=0; i<n; i++) temp_x[i] = SampleLW();
        
        for(int i=0; i<n; i++) v1[i] += (int128_t)temp_x[i] * b_pow;
        

        const int64_t *c_ptr = key->migd_key.c_coeffs + j*n;
        poly_mul_adj_acc_128(v2, c_ptr, temp_x, n);
        
        for(int i=0; i<n; i++) temp_x[i] = SampleLW();
        for(int i=0; i<n; i++) v2[i] += (int128_t)temp_x[i] * b_pow;

        b_pow *= b;
    }
    
    for(int r=0; r<MIGD_ROWS; r++) {
        get_eigd_poly_from_ctx(&key->migd_key.x_ctx, r, temp_poly, n);
        
        for(int i=0; i<n; i++) temp_x[i] = SampleLW();
        
        poly_mul_acc_128(v1, temp_poly, temp_x, n);
    }
    
    for(int r=0; r<MIGD_ROWS; r++) {
        get_eigd_poly_from_ctx(&key->migd_key.y_ctx, r, temp_poly, n);
        
        for(int i=0; i<n; i++) temp_x[i] = SampleLW();
        
        poly_mul_acc_128(v2, temp_poly, temp_x, n);
    }

#ifdef ZITAKA_DEBUG
    {
        char buf1[64], buf2[64];
        int128_t max_v1 = get_real_max_abs_128(v1, n);
        int128_t max_v2 = get_real_max_abs_128(v2, n);
        
        int128_to_str(max_v1, buf1);
        int128_to_str(max_v2, buf2);
        
        ZITAKA_LOG("[Sample] Offline Accumulators: Max(|v1|): %s, Max(|v2|): %s", buf1, buf2);
    }
#endif

    uint64_t den = 1ULL << 63; 
    
    for(int i=0; i<n; i++) {

        out_p1[i] = (int64_t)SampleArbitraryCenter128(v1[i], den);
        out_p2[i] = (int64_t)SampleArbitraryCenter128(v2[i], den);
    }

cleanup:
    free(v1); free(v2); 
    free(temp_x); free(temp_poly);
}


void OnlineSamp(const int64_t *u_hat, 
                const int128_t *c1_num, const int128_t *c2_num,
                int64_t *out_z1, int64_t *out_z2) {
    
    int n = ANTRAG_D; // 512
    
    uint64_t p = (1ULL << 28);
    uint64_t den = p * ANTRAG_Q; 

    for (int i = 0; i < n; i++) {
        out_z2[i] = (int64_t)SampleArbitraryCenter128(c2_num[i], den);
    }
#ifdef ZITAKA_DEBUG
    ZITAKA_LOG("[Sample] Online z2 sampled. z2[0]: %ld", out_z2[0]);
#endif

    for (int i = 0; i < n; i++) {
 
        int128_t conv_sum = 0;

        for (int j = 0; j <= i; j++) {
            conv_sum += (int128_t)out_z2[j] * u_hat[i - j];
        }

        for (int j = i + 1; j < n; j++) {
            conv_sum -= (int128_t)out_z2[j] * u_hat[i - j + n];
        }

        int128_t c1_prime_val = c1_num[i] - (conv_sum * ANTRAG_Q);
        out_z1[i] = (int64_t)SampleArbitraryCenter128(c1_prime_val, den);
#ifdef ZITAKA_DEBUG
        if (i == 0) {
            double center = (double)c1_prime_val / (double)den;
            double expected_part1 = (double)c1_num[i] / (double)den;
            double expected_part2 = (double)(conv_sum * ANTRAG_Q) / (double)den;
            ZITAKA_LOG("[OnlineSamp] i=0 Analysis:");
            ZITAKA_LOG("  c1 term (c1_num/den) = %f", expected_part1);
            ZITAKA_LOG("  u  term (conv/p)     = %f", expected_part2);
            ZITAKA_LOG("  Final Center         = %f", center);
        }
#endif
    }
}