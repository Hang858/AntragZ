/* src/sign.c */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h> 
#include "common.h"
#include "operator_interface.h"
#include "zalcon_samp.h"
#include "prematrix.h" 
#include "hash.h"
#include "rng.h"

// 必须使用 128 位整数来存储中间计算结果
typedef __int128_t int128_t;

// 声明 Compress 函数
size_t compress_sig(uint8_t *out, size_t max_out_len, const int16_t *x);

// 外部函数声明 
void OfflineSamp(const PreMatrix_Output *key, int64_t *out_p1, int64_t *out_p2);
void OnlineSamp(const int64_t *u_hat, 
                const int128_t *c1_num, const int128_t *c2_num, 
                int64_t *out_z1, int64_t *out_z2);

// ==========================================
// 调试辅助函数：打印统计信息
// ==========================================
static void print_stats_int64(const char* name, const int64_t* v, int n) {
    int64_t min_v = v[0], max_v = v[0];
    double sum = 0, sum_sq = 0;
    for(int i=0; i<n; i++) {
        if(v[i] < min_v) min_v = v[i];
        if(v[i] > max_v) max_v = v[i];
        sum += v[i];
        sum_sq += (double)v[i] * v[i];
    }
    printf("[DEBUG] %-10s: Min=%ld, Max=%ld, Avg=%.2f, RMS=%.2e\n", 
           name, min_v, max_v, sum/n, sqrt(sum_sq/n));
}

static void print_stats_int128(const char* name, const int128_t* v, int n) {
    // 为了打印方便，我们转换成 double (会丢失精度但足够看量级)
    double min_v = (double)v[0], max_v = (double)v[0];
    double sum_abs = 0;
    for(int i=0; i<n; i++) {
        double val = (double)v[i];
        if(val < min_v) min_v = val;
        if(val > max_v) max_v = val;
        sum_abs += fabs(val);
    }
    printf("[DEBUG] %-10s: Min=%.2e, Max=%.2e, AvgAbs=%.2e (128-bit)\n", 
           name, min_v, max_v, sum_abs/n);
}

// ==========================================
// 辅助函数
// ==========================================

// 1. (int64) += (int8) * (int64)
static void poly_mul_int8_int64_acc(int64_t *res, const int8_t *a, const int64_t *b, int n) {
    for (int i = 0; i < n; i++) {
        int8_t ai = a[i];
        if (ai == 0) continue;
        for (int j = 0; j < n; j++) {
            int64_t val = (int64_t)ai * b[j];
            int k = i + j;
            if (k < n) res[k] += val;
            else       res[k - n] -= val;
        }
    }
}

// 2. (int128) += (int8) * (int64)
static void poly_mul_int8_int64_to_128_acc(int128_t *res, const int8_t *a, const int64_t *b, int n) {
    for (int i = 0; i < n; i++) {
        int8_t ai = a[i];
        if (ai == 0) continue;
        for (int j = 0; j < n; j++) {
            int128_t val = (int128_t)ai * b[j];
            int k = i + j;
            if (k < n) res[k] += val;
            else       res[k - n] -= val;
        }
    }
}

// 3. (int128) += (int64) * (int64)
static void poly_mul_acc_128_64(int128_t *res, const int64_t *a, const int64_t *b, int n) {
    for (int i = 0; i < n; i++) {
        int64_t ai = a[i];
        if (ai == 0) continue;
        for (int j = 0; j < n; j++) {
            int128_t val = (int128_t)ai * b[j];
            int k = i + j;
            if (k < n) res[k] += val;
            else       res[k - n] -= val;
        }
    }
}

// ==========================================
// 4. 签名主函数 Algorithm 4
// ==========================================
int crypto_sign(uint8_t *sig, size_t *sig_len, const uint8_t *m, size_t mlen, const PrivateKey *sk) {
    int n = ANTRAG_D;
    
    // 1. Salt
    uint8_t salt[40];
    for (int i = 0; i < 40; i++) salt[i] = (uint8_t)(get_secure_random_u64() & 0xFF);

    // 2. Hash
    int16_t *c_hash = malloc(n * sizeof(int16_t));
    if (!c_hash) return -1;
    hash_to_point(c_hash, salt, 40, m, mlen);
    
    // Debug: Hash Stats
    {
        int64_t c_h_64[ANTRAG_D];
        for(int i=0;i<n;i++) c_h_64[i] = c_hash[i];
        print_stats_int64("c_hash", c_h_64, n);
    }

    // 3. Offline Sampling
    int64_t *p1 = calloc(n, sizeof(int64_t));
    int64_t *p2 = calloc(n, sizeof(int64_t));
    OfflineSamp(&sk->mat, p1, p2);
    
    print_stats_int64("p1 (off)", p1, n);
    print_stats_int64("p2 (off)", p2, n);

    // 4. Linear Transform
    int64_t *v1 = p1; 
    for (int i = 0; i < n; i++) v1[i] = -p1[i];

    int64_t *v2 = p2; 
    for (int i = 0; i < n; i++) v2[i] = (int64_t)c_hash[i] - p2[i]; // v2 = c - p2

    // Calc T = f*v2 - g*v1
    int64_t *T = calloc(n, sizeof(int64_t));
    poly_mul_int8_int64_acc(T, sk->f, v2, n);
    
    int8_t *neg_g = malloc(n * sizeof(int8_t));
    for (int i = 0; i < n; i++) neg_g[i] = -sk->g[i];
    poly_mul_int8_int64_acc(T, neg_g, v1, n);
    free(neg_g);
    
    print_stats_int64("T", T, n); // T 不应过大

    // c2_in = T * p
    int64_t p_val = 1LL << 28;
    int128_t *c2_in = malloc(n * sizeof(int128_t));
    for (int i = 0; i < n; i++) c2_in[i] = (int128_t)T[i] * p_val;

    // c1_in calculation
    int128_t *c1_in = calloc(n, sizeof(int128_t));
    
    poly_mul_int8_int64_to_128_acc(c1_in, sk->G, v1, n);
    
    int8_t *neg_F = malloc(n * sizeof(int8_t));
    for (int i = 0; i < n; i++) neg_F[i] = -sk->F[i];
    poly_mul_int8_int64_to_128_acc(c1_in, neg_F, v2, n);
    free(neg_F);

    for (int i = 0; i < n; i++) c1_in[i] *= p_val;

    // part2 = u_hat * T
    poly_mul_acc_128_64(c1_in, sk->mat.u_hat_num, T, n);
    
    print_stats_int128("c1_in", c1_in, n);
    print_stats_int128("c2_in", c2_in, n);

    free(T);
    free(p1); free(p2); 

    // 5. Online Sampling
    int64_t *z1 = malloc(n * sizeof(int64_t));
    int64_t *z2 = malloc(n * sizeof(int64_t));
    
    OnlineSamp(sk->mat.u_hat_num, c1_in, c2_in, z1, z2);
    
    print_stats_int64("z1", z1, n);
    print_stats_int64("z2", z2, n);
    
    free(c1_in); free(c2_in);

    // 6. Compute Signature s = c - B*z
    // s1 = -(f*z1 + F*z2)
    // s2 = c_hash - (g*z1 + G*z2)
    int16_t *s1 = malloc(n * sizeof(int16_t));
    int64_t *tmp_acc = calloc(n, sizeof(int64_t)); 

    // Calc s1
    poly_mul_int8_int64_acc(tmp_acc, sk->f, z1, n);
    poly_mul_int8_int64_acc(tmp_acc, sk->F, z2, n);
    
    int s1_overflow = 0;
    for (int i = 0; i < n; i++) {
        int64_t val = -tmp_acc[i];
        if (val < -32768 || val > 32767) s1_overflow++;
        s1[i] = (int16_t)val;
    }
    print_stats_int64("s1_raw", tmp_acc, n); // 注意这里打印的是 -(s1)

    // Calc s2
    memset(tmp_acc, 0, n * sizeof(int64_t));
    poly_mul_int8_int64_acc(tmp_acc, sk->g, z1, n);
    poly_mul_int8_int64_acc(tmp_acc, sk->G, z2, n);
    
    int64_t norm_sq = 0;
    for (int i = 0; i < n; i++) {
        int64_t val_s2 = (int64_t)c_hash[i] - tmp_acc[i];
        norm_sq += (int64_t)s1[i] * s1[i] + val_s2 * val_s2;
    }
    
    // 调试 s2 的部分值
    printf("[DEBUG] s2_sample[0] = %ld (c=%d, Bz=%ld)\n", 
           (int64_t)c_hash[0] - tmp_acc[0], c_hash[0], tmp_acc[0]);

    free(tmp_acc); free(z1); free(z2); free(c_hash);

    // 7. Norm Check
    printf("[DEBUG] NormSq: %ld (Limit: 34715664)\n", norm_sq);
    if (s1_overflow > 0) printf("[DEBUG] WARN: s1 had %d coeffs overflow int16\n", s1_overflow);

    if (norm_sq > 34715664LL) {
        free(s1);
        return 0; 
    }

    // 8. Compress
    if (*sig_len < 40 + 1) { free(s1); return -1; }
    memcpy(sig, salt, 40);
    size_t comp_len = compress_sig(sig + 40, *sig_len - 40, s1);
    
    free(s1);

    if (comp_len == 0) {
        printf("[DEBUG] Compression FAILED.\n");
        return 0; 
    }

    return 1; // 成功

err:
    if (c_hash) free(c_hash);
    return -1;
}