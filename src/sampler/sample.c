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
static int get_bit_width_128(int128_t val) {
    if (val == 0) return 0;
    if (val < 0) val = -val;
    int bits = 0;
    unsigned __int128 uval = (unsigned __int128)val;
    while (uval > 0) {
        bits++;
        uval >>= 1;
    }
    return bits;
}

static int64_t get_max_abs_64(const int64_t *v, int n) {
    int64_t max = 0;
    for(int i=0; i<n; i++) {
        int64_t val = v[i] < 0 ? -v[i] : v[i];
        if(val > max) max = val;
    }
    return max;
}

// 获取 int16 数组的最大绝对值
static int16_t get_max_abs_16(const int16_t *v, int n) {
    int16_t max = 0;
    for(int i=0; i<n; i++) {
        int16_t val = v[i] < 0 ? -v[i] : v[i];
        if(val > max) max = val;
    }
    return max;
}

// 计算大致位宽 (Bit Width)
static int get_bit_width_64(int64_t val) {
    if (val < 0) val = -val;
    if (val == 0) return 0;
    int bits = 0;
    while ((1ULL << bits) <= (uint64_t)val) {
        bits++;
    }
    return bits;
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

// #ifdef ZITAKA_DEBUG
//     {
//         char buf1[64], buf2[64];
//         int128_t max_v1 = get_real_max_abs_128(v1, n);
//         int128_t max_v2 = get_real_max_abs_128(v2, n);
        
//         int128_to_str(max_v1, buf1);
//         int128_to_str(max_v2, buf2);
        
//         ZITAKA_LOG("[Sample] Offline Accumulators: Max(|v1|): %s, Max(|v2|): %s", buf1, buf2);
//     }
// #endif

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
                    const int64_t *c1_num, const int64_t *c2_num,
                    int16_t *out_z1, int16_t *out_z2) {
    
    int n = ANTRAG_D; // 512
    int64_t p_shift = 28; // p = 2^28
    int64_t den = (1LL << p_shift) * ANTRAG_Q;

// #ifdef ZITAKA_DEBUG
//     // [调试 1] 分析输入参数 c1, c2 的大小
//     int64_t max_c1 = get_max_abs_64(c1_num, n);
//     int64_t max_c2 = get_max_abs_64(c2_num, n);
//     int bits_c1 = get_bit_width_64(max_c1);
//     int bits_c2 = get_bit_width_64(max_c2);
    
//     ZITAKA_LOG("========== [OnlineSamp Analysis Start] ==========");
//     ZITAKA_LOG("Input c1_num : Max(|v|) = %ld (approx %d bits)", max_c1, bits_c1);
//     ZITAKA_LOG("Input c2_num : Max(|v|) = %ld (approx %d bits)", max_c2, bits_c2);
// #endif

    // 1. Sample z2
    for (int i = 0; i < n; i++) {
        out_z2[i] = (int16_t)SampleArbitraryCenter64(c2_num[i], den);
    }
// #ifdef ZITAKA_DEBUG
//     // [调试 2] 分析采样结果 z2
//     int16_t max_z2 = get_max_abs_16(out_z2, n);
//     ZITAKA_LOG("Output z2    : Max(|v|) = %d", max_z2);
// #endif

//     // 为了统计 c1_prime 的最大值，我们需要一个临时数组
// #ifdef ZITAKA_DEBUG
//     int64_t *debug_c1_prime = malloc(n * sizeof(int64_t));
//     int64_t max_conv_term = 0; // 统计卷积项的最大值
// #endif

    // 2. Sample z1 (Convolution using int64)
    for (int i = 0; i < n; i++) {
        int64_t conv_sum = 0;

        for (int j = 0; j <= i; j++) {
            conv_sum += (int64_t)out_z2[j] * u_hat[i - j];
        }
        for (int j = i + 1; j < n; j++) {
            conv_sum -= (int64_t)out_z2[j] * u_hat[i - j + n];
        }
        int64_t c1_prime_val = c1_num[i] - (conv_sum * ANTRAG_Q);
// #ifdef ZITAKA_DEBUG
//         // 记录数据用于统计
//         int64_t term_to_sub = conv_sum * ANTRAG_Q;
//         debug_c1_prime[i] = c1_prime_val;
//         if (term_to_sub < 0) term_to_sub = -term_to_sub;
//         if (term_to_sub > max_conv_term) max_conv_term = term_to_sub;
// #endif
        out_z1[i] = (int16_t)SampleArbitraryCenter64(c1_prime_val, den);
    }
// #ifdef ZITAKA_DEBUG
//     // [调试 3] 分析中间变量 c1' 和卷积项，以及最终结果 z1
//     int64_t max_c1_prime = get_max_abs_64(debug_c1_prime, n);
//     int bits_c1_prime = get_bit_width_64(max_c1_prime);
//     int bits_conv = get_bit_width_64(max_conv_term);
//     int16_t max_z1 = get_max_abs_16(out_z1, n);

//     ZITAKA_LOG("Interm Conv*q: Max(|v|) = %ld (approx %d bits)", max_conv_term, bits_conv);
//     ZITAKA_LOG("Interm c1'   : Max(|v|) = %ld (approx %d bits)", max_c1_prime, bits_c1_prime);
//     ZITAKA_LOG("Output z1    : Max(|v|) = %d", max_z1);
    
//     // 安全性检查提示
//     if (bits_c1_prime > 62) {
//         ZITAKA_LOG("[WARNING] c1' bit width is dangerously close to 63 bits!");
//     } else {
//         ZITAKA_LOG("[PASS] Bit width safety check passed (<= 62 bits).");
//     }

//     ZITAKA_LOG("========== [OnlineSamp Analysis End] ==========");
    
//     free(debug_c1_prime);
// #endif
}