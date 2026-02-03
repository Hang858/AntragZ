#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "common.h"

// 声明 KeyGen 函数
extern int crypto_sign_keypair(PublicKey *pk, PrivateKey *sk);

// ==========================================
// 辅助工具：打印数组
// ==========================================
void print_int8_array(const char *label, const int8_t *data, int len, int limit) {
    printf("%s (前 %d 项): [", label, limit);
    for (int i = 0; i < limit && i < len; i++) {
        printf("%d", data[i]);
        if (i < limit - 1 && i < len - 1) printf(", ");
    }
    if (len > limit) printf(", ... (%d more)", len - limit);
    printf("]\n");
}

void print_int16_array(const char *label, const int16_t *data, int len, int limit) {
    printf("%s (前 %d 项): [", label, limit);
    for (int i = 0; i < limit && i < len; i++) {
        printf("%d", data[i]);
        if (i < limit - 1 && i < len - 1) printf(", ");
    }
    if (len > limit) printf(", ... (%d more)", len - limit);
    printf("]\n");
}

void print_int64_array(const char *label, const int64_t *data, int len, int limit) {
    printf("%s (前 %d 项): [", label, limit);
    for (int i = 0; i < limit && i < len; i++) {
        printf("%ld", data[i]);
        if (i < limit - 1 && i < len - 1) printf(", ");
    }
    if (len > limit) printf(", ... (%d more)", len - limit);
    printf("]\n");
}

// 辅助：取模
static inline int32_t mod_q(int32_t a) {
    int32_t r = a % ANTRAG_Q;
    return r < 0 ? r + ANTRAG_Q : r;
}

// ---------------------------------------------------------
// 测试 1: 验证公钥关系 h * f = g mod q
// ---------------------------------------------------------
int verify_public_key_relation(const PublicKey *pk, const PrivateKey *sk) {
    printf("\n[TEST] Verifying h * f == g mod q ...\n");
    
    int n = ANTRAG_D;
    int16_t *res = calloc(n, sizeof(int16_t));
    int16_t *f_poly = calloc(n, sizeof(int16_t));

    // f 转为 int16
    for(int i=0; i<n; i++) f_poly[i] = mod_q(sk->f[i]);

    // 计算 h * f
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int k = i + j;
            int32_t val = (int32_t)pk->h[i] * f_poly[j];
            if (k < n) 
                res[k] = mod_q(res[k] + val);
            else 
                res[k - n] = mod_q(res[k - n] - val);
        }
    }

    int fail = 0;
    for(int i=0; i<n; i++) {
        int16_t g_val = mod_q(sk->g[i]);
        if (res[i] != g_val) {
            if (fail < 5) printf("  Mismatch at %d: h*f=%d, g=%d\n", i, res[i], g_val);
            fail++;
        }
    }

    free(res);
    free(f_poly);

    if (fail) {
        printf("[FAIL] Public key relation check failed with %d errors.\n", fail);
        return 0;
    } else {
        printf("[PASS] Public key relation confirmed.\n");
        return 1;
    }
}

// ---------------------------------------------------------
// 测试 2: 验证 NTRU 方程
// ---------------------------------------------------------
int verify_ntru_equation(const PrivateKey *sk) {
    printf("\n[TEST] Verifying NTRU equation fG - gF = q ...\n");
    
    int n = ANTRAG_D;
    int64_t *res = calloc(n, sizeof(int64_t));

    // term1 = f * G
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            int k = i + j;
            int64_t val = (int64_t)sk->f[i] * sk->G[j];
            if(k < n) res[k] += val;
            else      res[k-n] -= val;
        }
    }

    // term2 = g * F
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            int k = i + j;
            int64_t val = (int64_t)sk->g[i] * sk->F[j];
            if(k < n) res[k] -= val;
            else      res[k-n] += val; 
        }
    }

    int fail = 0;
    int64_t q = ANTRAG_Q;

    if (res[0] != q) {
        // printf("  Mismatch at index 0: expected %ld, got %ld\n", q, res[0]);
        fail++;
    }

    for(int i=1; i<n; i++) {
        if (res[i] != 0) {
            // if(fail < 5) printf("  Mismatch at index %d: expected 0, got %ld\n", i, res[i]);
            fail++;
        }
    }

    free(res);

    if (fail) {
        printf("[FAIL] NTRU equation check failed.\n");
        return 0;
    } else {
        printf("[PASS] NTRU equation fG - gF = q confirmed.\n");
        return 1;
    }
}

// ---------------------------------------------------------
// 测试 3: 检查 PreMatrix 有效性
// ---------------------------------------------------------
int verify_prematrix(const PrivateKey *sk) {
    printf("\n[TEST] Verifying PreMatrix validity...\n");
    
    if (sk->mat.valid != 1) {
        printf("[FAIL] sk->mat.valid flag is 0.\n");
        return 0;
    }

    int nonzero = 0;
    for(int i=0; i<ANTRAG_D; i++) {
        if (sk->mat.c11[i] != 0) nonzero = 1;
    }

    if (!nonzero) {
        printf("[WARN] PreMatrix c11 is all zeros (suspicious).\n");
        return 0;
    }

    printf("[PASS] PreMatrix is flagged valid and contains data.\n");
    return 1;
}

// ==========================================
// 主函数
// ==========================================
int main() {
    printf("==========================================\n");
    printf("      Zitaka KeyGen Verification Tool     \n");
    printf("==========================================\n");

    PublicKey pk;
    PrivateKey sk;

    clock_t start = clock();
    int ret = crypto_sign_keypair(&pk, &sk);
    clock_t end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    if (!ret) {
        printf("[FATAL] crypto_sign_keypair returned 0 (Failure).\n");
        return 1;
    }

    printf("[INFO] KeyGen completed in %.4f seconds.\n", cpu_time_used);

    // === 新增：打印密钥信息 ===
    printf("\n=== Generated Keys Dump ===\n");
    
    // 1. 打印公钥 h
    print_int16_array("Public Key h", pk.h, ANTRAG_D, 32);

    // 2. 打印私钥 f, g, F, G
    print_int8_array("Secret Key f", sk.f, ANTRAG_D, 32);
    print_int8_array("Secret Key g", sk.g, ANTRAG_D, 32);
    print_int8_array("Secret Key F", sk.F, ANTRAG_D, 32);
    print_int8_array("Secret Key G", sk.G, ANTRAG_D, 32);

    // 3. 打印部分预矩阵信息 (验证是否有数据)
    print_int64_array("PreMatrix c11", sk.mat.c11, ANTRAG_D, 16);
    print_int64_array("PreMatrix c21", sk.mat.c21, ANTRAG_D, 16);
    
    printf("===========================\n");

    // 开始验证
    int score = 0;
    score += verify_public_key_relation(&pk, &sk);
    score += verify_ntru_equation(&sk);
    score += verify_prematrix(&sk);

    printf("\n------------------------------------------\n");
    if (score == 3) {
        printf("RESULT: ALL TESTS PASSED. KeyGen is correct.\n");
        return 0;
    } else {
        printf("RESULT: %d/3 TESTS PASSED. KeyGen has errors.\n", score);
        return 1;
    }
}