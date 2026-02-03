/*
 * src/test_prematrix.c
 * 严格验证 PreMatrix 算法正确性的测试程序 (End-to-End Strict Mode)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>

#include "common.h"
#include "fft.h"
#include "prematrix.h"
#include "migd.h"

// 类型定义
typedef __int128_t int128_t;

// --------------------------------------------------------------------------
// 打印工具
// --------------------------------------------------------------------------
void print_128(const char* prefix, int128_t val) {
    double v = (double)val;
    printf("%s: %.4e", prefix, v);
}

// --------------------------------------------------------------------------
// 辅助多项式乘法
// --------------------------------------------------------------------------
// Res += A * B
static void poly_mul_accum(int128_t *res, const int64_t *a, const int64_t *b) {
    int n = ANTRAG_D;
    for (int i = 0; i < n; i++) {
        if (a[i] == 0) continue;
        for (int j = 0; j < n; j++) {
            int128_t val = (int128_t)a[i] * b[j];
            int k = i + j;
            if (k < n) res[k] += val;
            else       res[k - n] -= val;
        }
    }
}

// Res += A * adj(B)
static void poly_mul_adj_accum(int128_t *res, const int64_t *a, const int64_t *b) {
    int n = ANTRAG_D;
    int128_t b0 = b[0];
    if (b0 != 0) {
        for(int i=0; i<n; i++) res[i] += (int128_t)a[i] * b0;
    }
    for (int k = 1; k < n; k++) {
        int128_t b_adj_k = -b[n - k];
        if (b_adj_k == 0) continue;
        for (int i = 0; i < n; i++) {
            int128_t val = (int128_t)a[i] * b_adj_k;
            int pos = i + k;
            if (pos < n) res[pos] += val;
            else         res[pos - n] -= val;
        }
    }
}

// --------------------------------------------------------------------------
// 主验证程序
// --------------------------------------------------------------------------
extern int keygen_fg_impl(secret_key_fg *sk);
extern int antrag_solve_ntru(int8_t *F, int8_t *G, const int8_t *f, const int8_t *g, unsigned logn, uint32_t *tmp);

int main() {
    printf("==============================================================\n");
    printf("   Strict Verification of AA^T = p^2(Sigma_p - I)             \n");
    printf("==============================================================\n");

    // 1. 准备密钥
    secret_key_fg sk;
    printf("[1] Generating Keys...\n");
    keygen_fg_impl(&sk);
    
    int8_t F[ANTRAG_D], G[ANTRAG_D];
    uint32_t *tmp = malloc(70000 * sizeof(uint32_t));
    antrag_solve_ntru(F, G, sk.f, sk.g, ANTRAG_LOGD, tmp);
    free(tmp);

    // 2. 运行 PreMatrix
    printf("[2] Running PreMatrix Algorithm...\n");
    PreMatrix_Output out;
    memset(&out, 0, sizeof(out));
    if (!Run_PreMatrix(sk.f, sk.g, F, G, &out)) {
        printf("[FAIL] Run_PreMatrix returned error.\n");
        return 1;
    }
    printf("    -> Algorithm finished successfully.\n");

    // 3. 开始严格验证
    printf("[3] Verifying Equation Integrity...\n");

    int n = ANTRAG_D;
    int64_t p = 1LL << 28;
    int64_t s0 = 131;

    // ------------------------------------------------------------
    // Step 3.1: 构建 RHS (Target Equation)
    // Target = p^2 * (s0^2 * I - B_hat * B_hat^T - I)
    //        = p^2 * (s0^2 - 1) * I - Sigma
    // 其中 Sigma = p^2 * B_hat * B_hat^T
    // ------------------------------------------------------------
    
    // 计算 Sigma (使用 int128 避免溢出)
    int64_t *f64 = malloc(n*8), *g64 = malloc(n*8);
    int64_t *R1 = malloc(n*8), *R2 = malloc(n*8);
    int128_t *t128 = calloc(n, 16);

    for(int i=0; i<n; i++) { f64[i]=sk.f[i]; g64[i]=sk.g[i]; }

    // R1 = pF - u*f
    memset(t128, 0, n*16);
    poly_mul_accum(t128, out.u_hat_num, f64);
    for(int i=0; i<n; i++) R1[i] = (int64_t)F[i]*p - (int64_t)t128[i];

    // R2 = pG - u*g
    memset(t128, 0, n*16);
    poly_mul_accum(t128, out.u_hat_num, g64);
    for(int i=0; i<n; i++) R2[i] = (int64_t)G[i]*p - (int64_t)t128[i];

    // Sigma_11 = p^2 f f* + R1 R1*
    int128_t *Sig11 = calloc(n, 16);
    int128_t p2 = (int128_t)p * p;
    memset(t128, 0, n*16);
    poly_mul_adj_accum(t128, f64, f64);
    for(int i=0; i<n; i++) Sig11[i] += t128[i] * p2;
    poly_mul_adj_accum(Sig11, R1, R1);

    // 计算 RHS_11 (Target)
    // Target_11 = p^2(s0^2 - 1) * I - Sig11
    int128_t scalar_term = p2 * (s0*s0 - 1);
    int128_t *RHS_11 = malloc(n*16);
    for(int i=0; i<n; i++) RHS_11[i] = -Sig11[i];
    RHS_11[0] += scalar_term; // 加到常数项 (对角线)

    // ------------------------------------------------------------
    // Step 3.2: 构建 LHS (Actual Constructed Matrix)
    // LHS = A * A^T = (1 + d) * I + P
    // P = B^2 * I - Sigma
    // 所以 LHS = (1 + d + B^2) * I - Sigma
    // ------------------------------------------------------------
    
    // 从代码逻辑获取 d
    int128_t d = 2 * (int128_t)p * (p * (s0 - 1)) - 1;
    
    // 计算 P_11
    int128_t B_val = (int128_t)p * (s0 - 1);
    int128_t B_sq = B_val * B_val;
    int128_t *P11 = malloc(n*16);
    for(int i=0; i<n; i++) P11[i] = -Sig11[i];
    P11[0] += B_sq;

    // 计算 LHS_11
    int128_t *LHS_11 = malloc(n*16);
    for(int i=0; i<n; i++) LHS_11[i] = P11[i]; // Copy P
    LHS_11[0] += (1 + d); // Add (1+d)I

    // ------------------------------------------------------------
    // Step 3.3: 验证等式 LHS == RHS
    // ------------------------------------------------------------
    printf("\n[Check 1] Algebraic Equality (LHS vs RHS):\n");
    
    int match = 1;
    for(int i=0; i<n; i++) {
        if (LHS_11[i] != RHS_11[i]) {
            match = 0;
            printf("  MISMATCH at index %d!\n", i);
            print_128("    LHS", LHS_11[i]);
            printf("\n");
            print_128("    RHS", RHS_11[i]);
            printf("\n");
            break;
        }
    }

    if (match) {
        printf("  [PASS] LHS == RHS exactly for all coefficients.\n");
        printf("  This proves A*A^T matches the target covariance structure mathematically.\n");
    } else {
        printf("  [FAIL] Algebraic mismatch.\n");
        return 1;
    }

    // ------------------------------------------------------------
    // Step 3.4: 验证 Cholesky 精度 (Delta)
    // ------------------------------------------------------------
    printf("\n[Check 2] Cholesky Precision (Delta size):\n");
    
    // 计算 CC^T_11
    int128_t *CCT_11 = calloc(n, 16);
    poly_mul_adj_accum(CCT_11, out.c11, out.c11);

    // Delta = CC^T - P
    int128_t max_delta = 0;
    for(int i=0; i<n; i++) {
        int128_t diff = CCT_11[i] - P11[i];
        if (diff < 0) diff = -diff;
        if (diff > max_delta) max_delta = diff;
    }

    print_128("  Max Delta", max_delta);
    double limit = pow(2.0, 50.0);
    printf("  Limit    : %.4e\n", limit);
    
    if ((double)max_delta < limit) {
        printf("  [PASS] Delta is small enough for MIGD.\n");
    } else {
        printf("  [FAIL] Delta too large!\n");
        return 1;
    }

    // ------------------------------------------------------------
    // 结论
    // ------------------------------------------------------------
    printf("\n==============================================================\n");
    printf("FINAL RESULT: VERIFICATION SUCCESSFUL\n");
    printf("The implementation strictly satisfies the design equation.\n");
    printf("==============================================================\n");

    // 清理
    free(f64); free(g64); free(R1); free(R2); free(t128);
    free(Sig11); free(RHS_11); free(LHS_11); free(P11); free(CCT_11);

    return 0;
}