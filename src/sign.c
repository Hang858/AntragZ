#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h> 
#include "common.h"
#include "operator_interface.h"
#include "zalcon_samp.h"
#include "sample.h"
#include "prematrix.h" 
#include "hash.h"
#include "rng.h"
#include "poly.h"
#include "encode.h"
#include "bench.h"

#ifdef ZITAKA_DEBUG
static int32_t get_max_abs_32(const int32_t *v, int n) {
    int32_t max = 0;
    for(int i=0; i<n; i++) {
        int32_t val = v[i] < 0 ? -v[i] : v[i];
        if(val > max) max = val;
    }
    return max;
}

static int16_t get_max_abs_16(const int16_t *v, int n) {
    int16_t max = 0;
    for(int i=0; i<n; i++) {
        int16_t val = v[i] < 0 ? -v[i] : v[i];
        if(val > max) max = val;
    }
    return max;
}
#endif

int crypto_sign(uint8_t *sig, size_t *sig_len, const uint8_t *m, size_t mlen, const PrivateKey *sk) {

    BENCH_INIT(total);
    BENCH_INIT(sample_off);
    BENCH_INIT(hash);
    BENCH_INIT(linear);
    BENCH_INIT(sample_on);
    BENCH_INIT(calc_s);
    BENCH_INIT(compress);

    BENCH_START(total);
    int n = ANTRAG_D;

    // 生成salt
    uint8_t salt[40];
    for (int i = 0; i < 40; i++) salt[i] = (uint8_t)(get_secure_random_u64() & 0xFF);
    BENCH_START(sample_off);
    // 离线采样生成 p1, p2
    int64_t *p1 = calloc(n, sizeof(int64_t));
    int64_t *p2 = calloc(n, sizeof(int64_t));
    OfflineSamp(&sk->mat, p1, p2);
    int16_t *v1 = malloc(n * sizeof(int16_t));
    int16_t *v2 = malloc(n * sizeof(int16_t));
    for(int i=0; i<n; i++) v1[i] = (int16_t)(-p1[i]); // 直接存 -p1
    for(int i=0; i<n; i++) v2[i] = (int16_t)(-p2[i]); // 暂存 -p2
    free(p1); 
    free(p2);
    BENCH_STOP(sample_off);



    BENCH_START(hash);
    // 在线消息哈希到 c 中
    int16_t *c_hash = malloc(n * sizeof(int16_t));
    if (!c_hash) return -1;
    hash_to_point(c_hash, salt, 40, m, mlen);
    BENCH_STOP(hash);

    ZITAKA_LOG("Hash point generated. c_hash[0]: %d", c_hash[0]);
    // v1, v2 构成了 c-p
    BENCH_START(linear);
    // 计算 f * v2 组成B_hat ^(-1) * (c - p)
    for (int i = 0; i < n; i++) v2[i] = c_hash[i] + v2[i];
    int32_t *T = calloc(n, sizeof(int32_t));
    poly_mul_int8_int16_to_32_acc(T, sk->f, v2, n);
    // 计算 -g * v1
    int8_t *neg_g = malloc(n * sizeof(int8_t));
    for (int i = 0; i < n; i++) neg_g[i] = -sk->g[i];
    poly_mul_int8_int16_to_32_acc(T, neg_g, v1, n);  // 此时 T 中为 f * v2 - g * v1
    free(neg_g);

    // c2_in = T * p 
    // (因为算法输入 u_hat 只存了整数部分，实际的 u_hat 为 存储的 u_hat 除以 p，为了让其他部分是同一个数量级，也乘p)
    int64_t p_val = 1LL << 28;
    int64_t *c2_in = malloc(n * sizeof(int64_t));
    for (int i = 0; i < n; i++) c2_in[i] = (int64_t)T[i] * p_val;

    // c1_in = u_hat (f * v2 - g * v1) + G * v1 - F * v2 = u_hat * T + (G * v1) - (F * v2)
    int64_t *c1_in = calloc(n, sizeof(int64_t));
    int32_t *tmp_small_acc = calloc(n, sizeof(int32_t));
    poly_mul_int8_int16_to_32_acc(tmp_small_acc, sk->G, v1, n);
    
    int8_t *neg_F = malloc(n * sizeof(int8_t));
    for (int i = 0; i < n; i++) neg_F[i] = -sk->F[i];
    poly_mul_int8_int16_to_32_acc(tmp_small_acc, neg_F, v2, n);
    free(neg_F);

    for(int i=0; i<n; i++) c1_in[i] = (int64_t)tmp_small_acc[i] * p_val;
    free(tmp_small_acc);
    // part2 = u_hat * T
    poly_mul_int64_int32_to_64_acc(c1_in, sk->mat.u_hat_num, T, n);

    free(T);
    free(v1); free(v2);
    BENCH_STOP(linear);
    ZITAKA_LOG("Starting OnlineSamp...");
    BENCH_START(sample_on);
    // Online Sampling
    int16_t *z1 = malloc(n * sizeof(int16_t));
    int16_t *z2 = malloc(n * sizeof(int16_t));
    
    OnlineSamp(sk->mat.u_hat_num, c1_in, c2_in, z1, z2);
    free(c1_in); free(c2_in);
    BENCH_STOP(sample_on);
    BENCH_START(calc_s);
    int32_t *final_acc = calloc(n, sizeof(int32_t)); 
    // 计算 c - (B * Z)
    poly_mul_int8_int16_to_32_acc(final_acc, sk->f, z1, n);
    poly_mul_int8_int16_to_32_acc(final_acc, sk->F, z2, n);
// #ifdef ZITAKA_DEBUG
//     int32_t max_Bz1 = get_max_abs_32(final_acc, n);
// #endif

    int16_t *s1 = malloc(n * sizeof(int16_t));
    for (int i = 0; i < n; i++) s1[i] = (int16_t)(-final_acc[i]); // s1 = -(B*z)_1

    memset(final_acc, 0, n * sizeof(int32_t));
    poly_mul_int8_int16_to_32_acc(final_acc, sk->g, z1, n);
    poly_mul_int8_int16_to_32_acc(final_acc, sk->G, z2, n);
// #ifdef ZITAKA_DEBUG
//     int32_t max_Bz2 = get_max_abs_32(final_acc, n);
//     int16_t max_c = get_max_abs_16(c_hash, n);
// #endif

    free(z1); free(z2);

// #ifdef ZITAKA_DEBUG
//     int32_t max_s2 = 0;
// #endif
    
    int64_t norm_sq = 0;
    for (int i = 0; i < n; i++) {
        // s2 = c - (g*z1 + G*z2)
        int32_t val_s2 = (int32_t)c_hash[i] - final_acc[i];
        norm_sq += (int64_t)s1[i] * s1[i] + (int64_t)val_s2 * val_s2;
    }
// #ifdef ZITAKA_DEBUG
//     ZITAKA_LOG("========== [Magnitude Analysis] ==========");
//     ZITAKA_LOG("1. Target c (hash): Max(|c2|) = %d", max_c);
//     ZITAKA_LOG("2. Vector B*z     : Max(|Bz1|) = %d, Max(|Bz2|) = %d", max_Bz1, max_Bz2);
//     ZITAKA_LOG("   -> Note: Bz2 should be very close to c2");
//     ZITAKA_LOG("3. Signature s    : Max(|s1|) = %d, Max(|s2|) = %d", get_max_abs_16(s1, n), max_s2);
//     ZITAKA_LOG("4. Norm Squared   : %ld (Limit: 34715664)", norm_sq);
//     ZITAKA_LOG("==========================================");
// #endif

    free(final_acc);
    free(c_hash);
    BENCH_STOP(calc_s);
    if (norm_sq > 34715664LL) {
        ZITAKA_LOG("REJECTED: Norm too large. Val: %ld > Limit: %ld", norm_sq, 34715664LL);
        free(s1);
        return 0; 
    }
    BENCH_START(compress);
    if (*sig_len < 40 + 1) { free(s1); return -1; }
    memcpy(sig, salt, 40);
    size_t comp_len = compress_sig(sig + 40, *sig_len - 40, s1);
    
    free(s1);
    BENCH_STOP(compress);

    BENCH_STOP(total);
    BENCH_PRINT(sample_off, "Offline Sampling");
    BENCH_PRINT(hash,       "Hashing");
    BENCH_PRINT(linear,     "Linear Transform");
    BENCH_PRINT(sample_on,  "Online Sampling");
    BENCH_PRINT(calc_s,     "Calc s = c - Bz");
    BENCH_PRINT(compress,   "Compression");
    BENCH_PRINT(total,      "Total Sign");
    if (comp_len == 0) {
        ZITAKA_LOG("Compression algorithm failed (buffer overflow?)");
        return 0; 
    }
    *sig_len = 40 + comp_len;
    ZITAKA_LOG("Sign Success. Final Size: %zu bytes", *sig_len);
    return 1;

err:
    if (c_hash) free(c_hash);
    return -1;
}