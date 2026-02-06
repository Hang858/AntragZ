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
static int64_t get_max_abs(const int64_t *v, int n) {
    int64_t max = 0;
    for(int i=0; i<n; i++) {
        int64_t val = v[i] < 0 ? -v[i] : v[i];
        if(val > max) max = val;
    }
    return max;
}

static int get_bit_width(int64_t val) {
    if (val < 0) val = -val;
    if (val == 0) return 0;
    int bits = 0;
    while ((1ULL << bits) <= (uint64_t)val) {
        bits++;
    }
    return bits;
}
#endif

int crypto_sign(uint8_t *sig, size_t *sig_len, const uint8_t *m, size_t mlen, const PrivateKey *sk) {

    BENCH_INIT(total);
    BENCH_INIT(sample_off);
    BENCH_INIT(hash);
    BENCH_INIT(linear);
    BENCH_INIT(sample_on);
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
    BENCH_STOP(sample_off);

#ifdef ZITAKA_DEBUG
    int64_t max_p1 = get_max_abs(p1, n);
    int64_t max_p2 = get_max_abs(p2, n);
    int bits_p1 = get_bit_width(max_p1);
    int bits_p2 = get_bit_width(max_p2);
    
    ZITAKA_LOG("[Stats] Offline p-vectors:");
    ZITAKA_LOG("  Max(|p1|): %ld (approx %d bits)", max_p1, bits_p1);
    ZITAKA_LOG("  Max(|p2|): %ld (approx %d bits)", max_p2, bits_p2);
    
    if (bits_p1 > 12 || bits_p2 > 12) {
        ZITAKA_LOG("  [WARNING] p-vector values are unusually large! Check sigma scaling.");
    }
#endif

#ifdef ZITAKA_DEBUG
    ZITAKA_LOG("Offline Sample done. Max(|p1|): %ld", get_max_abs(p1, n));
#endif
    BENCH_START(hash);
    // 在线消息哈希到 c 中
    int16_t *c_hash = malloc(n * sizeof(int16_t));
    if (!c_hash) return -1;
    hash_to_point(c_hash, salt, 40, m, mlen);
    BENCH_STOP(hash);

    ZITAKA_LOG("Hash point generated. c_hash[0]: %d", c_hash[0]);
    BENCH_START(linear);
    // 计算 -p1
    int64_t *v1 = p1; 
    for (int i = 0; i < n; i++) v1[i] = -p1[i];
    // 计算 H-p2
    int64_t *v2 = p2; 
    for (int i = 0; i < n; i++) v2[i] = (int64_t)c_hash[i] - p2[i]; // v2 = H - p2
    // v1, v2 构成了 c-p

    // 计算 f * v2 组成B_hat ^(-1) * (c - p)
    int64_t *T = calloc(n, sizeof(int64_t));
    poly_mul_int8_int64_acc(T, sk->f, v2, n);
    // 计算 -g * v1
    int8_t *neg_g = malloc(n * sizeof(int8_t));
    for (int i = 0; i < n; i++) neg_g[i] = -sk->g[i];
    poly_mul_int8_int64_acc(T, neg_g, v1, n);  // 此时 T 中为 f * v2 - g * v1
    free(neg_g);
    BENCH_STOP(linear);
#ifdef ZITAKA_DEBUG
    ZITAKA_LOG("Intermediate T computed. Max(|T|): %ld", get_max_abs(T, n));
#endif

    // c2_in = T * p 
    // (因为算法输入 u_hat 只存了整数部分，实际的 u_hat 为 存储的 u_hat 除以 p，为了让其他部分是同一个数量级，也乘p)
    int64_t p_val = 1LL << 28;
    int128_t *c2_in = malloc(n * sizeof(int128_t));
    for (int i = 0; i < n; i++) c2_in[i] = (int128_t)T[i] * p_val;

    // c1_in = u_hat (f * v2 - g * v1) + G * v1 - F * v2 = u_hat * T + (G * v1) - (F * v2)
    int128_t *c1_in = calloc(n, sizeof(int128_t));
    
    poly_mul_int8_int64_to_128_acc(c1_in, sk->G, v1, n);
    
    int8_t *neg_F = malloc(n * sizeof(int8_t));
    for (int i = 0; i < n; i++) neg_F[i] = -sk->F[i];
    poly_mul_int8_int64_to_128_acc(c1_in, neg_F, v2, n);
    free(neg_F);

    for (int i = 0; i < n; i++) c1_in[i] *= p_val;

    // part2 = u_hat * T
    poly_mul_acc_128(c1_in, sk->mat.u_hat_num, T, n);

    free(T);
    free(p1); free(p2); 

    ZITAKA_LOG("Starting OnlineSamp...");
    BENCH_START(sample_on);
    // Online Sampling
    int64_t *z1 = malloc(n * sizeof(int64_t));
    int64_t *z2 = malloc(n * sizeof(int64_t));
    
    OnlineSamp(sk->mat.u_hat_num, c1_in, c2_in, z1, z2);
    free(c1_in); free(c2_in);
    BENCH_STOP(sample_on);
#ifdef ZITAKA_DEBUG
    ZITAKA_LOG("Sampled z vectors. z1[0]: %ld, z2[0]: %ld", z1[0], z2[0]);
#endif

    int16_t *s1 = malloc(n * sizeof(int16_t));
    int64_t *tmp_acc = calloc(n, sizeof(int64_t)); 
    // 计算 c - (B * Z)
    poly_mul_int8_int64_acc(tmp_acc, sk->f, z1, n);
    poly_mul_int8_int64_acc(tmp_acc, sk->F, z2, n);

    int s1_overflow = 0;
    for (int i = 0; i < n; i++) {
        int64_t val = -tmp_acc[i];
        if (val < -32768 || val > 32767) s1_overflow++;
        s1[i] = (int16_t)val;
    }

    if (s1_overflow > 0) {
        ZITAKA_LOG("WARNING: s1 overflow detected count: %d (Bad sign)", s1_overflow);
    }

    memset(tmp_acc, 0, n * sizeof(int64_t));
    poly_mul_int8_int64_acc(tmp_acc, sk->g, z1, n);
    poly_mul_int8_int64_acc(tmp_acc, sk->G, z2, n);
    
    int64_t norm_sq = 0;
    for (int i = 0; i < n; i++) {
        int64_t val_s2 = (int64_t)c_hash[i] - tmp_acc[i];
        norm_sq += (int64_t)s1[i] * s1[i] + val_s2 * val_s2;
    }
    free(tmp_acc); free(z1); free(z2); free(c_hash);

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
    printf("=== Sign Performance (Cycles) ===\n");
    BENCH_PRINT(sample_off, "Offline Sampling");
    BENCH_PRINT(hash,       "Hashing");
    BENCH_PRINT(linear,     "Linear Transform");
    BENCH_PRINT(sample_on,  "Online Sampling");
    BENCH_PRINT(compress,   "Compression");
    BENCH_PRINT(total,      "Total Sign");
    printf("=================================\n");
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