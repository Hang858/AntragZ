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

static void get_eigd_poly_from_ctx(const MemContext *ctx, int row, int64_t *out_poly, int n) {
    for(int i=0; i<n; i++) {
        out_poly[i] = matrix_get(ctx, row, i);
    }
}

void OfflineSamp(const PreMatrix_Output *key, int64_t *out_p1, int64_t *out_p2) {
    int n = ANTRAG_D;
    
    int128_t *v1 = calloc(n, sizeof(int128_t));
    int128_t *v2 = calloc(n, sizeof(int128_t));
    
    int64_t *temp_x = malloc(n * sizeof(int64_t));
    int64_t *temp_poly = malloc(n * sizeof(int64_t));
    
    if (!v1 || !v2 || !temp_x || !temp_poly) goto cleanup;

    // =================================================
    // Block 1: Identity Matrix I (2 columns)
    // =================================================
    // Col 1: [1; 0] * x1 -> v1 += x1
    for(int i=0; i<n; i++) v1[i] += SampleLW(); 
    
    // Col 2: [0; 1] * x2 -> v2 += x2
    for(int i=0; i<n; i++) v2[i] += SampleLW();

    // =================================================
    // Block 2: Cholesky Matrix C = [c11 0; c21 c22]
    // =================================================
    // Col 1: [c11; c21] * x
    for(int i=0; i<n; i++) temp_x[i] = SampleLW();
    poly_mul_acc_128(v1, key->c11, temp_x, n); // v1 += c11 * x
    poly_mul_acc_128(v2, key->c21, temp_x, n); // v2 += c21 * x
    

    for(int i=0; i<n; i++) temp_x[i] = SampleLW();
    poly_mul_acc_128(v2, key->c22, temp_x, n); // v2 += c22 * x


    int64_t b_pow = 1;
    const int64_t b = 32768;
    
    for(int j=0; j<K_VAL; j++) {

        for(int i=0; i<n; i++) temp_x[i] = SampleLW();
        
        for(int i=0; i<n; i++) v1[i] += (int128_t)temp_x[i] * b_pow;
        

        const int64_t *c_ptr = key->migd_key.c_coeffs + j*n;
        poly_mul_adj_acc_128(v2, c_ptr, temp_x, n);
        
        // --- L_j Col 2: [0; b^j] * x ---
        for(int i=0; i<n; i++) temp_x[i] = SampleLW();
        
        // v2 += b^j * x
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
        // out_z2[i] = SampleSW(c2_num[i], den);
        // out_z2[i] = (int64_t)llround((double)c2_num[i] / (double)den);
        out_z2[i] = (int64_t)SampleArbitraryCenter128(c2_num[i], den);
        // if (i == 0) {
        //     printf("[DEBUG-Z2] i=0: center=%.2f, sampled=%ld, diff=%.2f\n", 
        //            (double)c2_num[i]/den, out_z2[i], (double)out_z2[i] - (double)c2_num[i]/den);
        // }
    }

    for (int i = 0; i < n; i++) {
 
        int128_t conv_sum = 0;

        for (int j = 0; j <= i; j++) {
            conv_sum += (int128_t)out_z2[j] * u_hat[i - j];
        }

        for (int j = i + 1; j < n; j++) {
            conv_sum -= (int128_t)out_z2[j] * u_hat[i - j + n];
        }

        int128_t c1_prime_val = c1_num[i] - (conv_sum * ANTRAG_Q);
        // out_z1[i] = SampleSW(c1_prime_val, den);
        // out_z1[i] = (int64_t)llround((double)c1_prime_val / (double)den);
        out_z1[i] = (int64_t)SampleArbitraryCenter128(c1_prime_val, den);
        // if (i == 0) {
        //     printf("[DEBUG-Z1] i=0: center_c1_prime=%.2f, sampled=%ld, diff=%.2f\n", 
        //            (double)c1_prime_val/den, out_z1[i], (double)out_z1[i] - (double)c1_prime_val/den);
        // }
        double total_diff_sq = 0;
        for(int i=0; i<n; i++) {
            double diff = (double)out_z1[i] - (double)c1_prime_val/den; // 这里需要逻辑上拿到当前的 c1_prime_val
            total_diff_sq += diff * diff;
        }
    }
}