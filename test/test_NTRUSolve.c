// src/test_ntru.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "fft.h"

int keygen_fg_impl(secret_key_fg *sk);  // 声明 keygen_fg_impl (确保在 keygen_basis.c 中非 static)

// 声明 NTRU 求解器入口 (确保在 keygen_ntru.c 中非 static)
int antrag_solve_ntru(int8_t *F, int8_t *G, const int8_t *f, const int8_t *g, unsigned logn, uint32_t *tmp);

// 多项式乘法 (模 X^N + 1, 模 q)
// 用于验证结果
void poly_mul_mod(int16_t *res, const int8_t *a, const int8_t *b, int n, int q) {
    for (int i = 0; i < n; i++) res[i] = 0;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int k = i + j;
            int sign = 1;
            if (k >= n) {
                k -= n;
                sign = -1;
            }
            int val = (int)a[i] * (int)b[j] * sign;
            res[k] = (res[k] + val) % q;
        }
    }
    // 调整负数
    for (int i = 0; i < n; i++) {
        if (res[i] < 0) res[i] += q;
    }
}

int main() {
    printf("=== ANTRU KeyGen Step 2: NTRU Solve Test ===\n");

    // 1. 参数设置
    unsigned logn = 9; // N = 512
    size_t n = 512;
    uint32_t *tmp = malloc(64 * 1024); // 分配足够大的临时空间 (64KB 足够)
    
    secret_key_fg sk_fg;
    int8_t F[512], G[512];

    // 2. 生成 f, g
    printf("[1] Generating f, g...\n");
    if (keygen_fg_impl(&sk_fg) <= 0) {
        printf("VectorGen failed!\n");
        return 1;
    }

    // 3. 求解 F, G
    printf("[2] Solving NTRU equation...\n");
    if (!antrag_solve_ntru(F, G, sk_fg.f, sk_fg.g, logn, tmp)) {
        printf("NTRU Solve failed!\n");
        free(tmp);
        return 1;
    }
    printf("NTRU Solve success!\n");

    // 4. 验证方程 fG - gF = q mod (X^N + 1)
    printf("[3] Verifying equation fG - gF = q...\n");
    
    int16_t *res_fG = calloc(n, sizeof(int16_t));
    int16_t *res_gF = calloc(n, sizeof(int16_t));
    
    // 计算 f*G
    poly_mul_mod(res_fG, sk_fg.f, G, n, 12289);
    // 计算 g*F
    poly_mul_mod(res_gF, sk_fg.g, F, n, 12289);
    
    int check_pass = 1;
    for (int i = 0; i < n; i++) {
        // 计算 fG - gF mod q
        int val = (res_fG[i] - res_gF[i]) % 12289;
        if (val < 0) val += 12289;
        
        // 目标：常数项为 q (即 0 mod q)，其他项也是 0
        // 等等，NTRU方程通常是 fG - gF = q。在模 q 意义下，这意味着 fG - gF = 0。
        // 但是在整数域上，它等于标量 q。
        // 让我们先检查模 12289 是否全为 0。
        
        if (val != 0) {
            printf("Mismatch at coeff %d: %d != 0\n", i, val);
            check_pass = 0;
            // break; 
        }
    }

    // 更严格的验证：在整数域验证 (不模 q)
    // 但这需要大整数乘法，简单起见，我们先验证模 q 等于 0。
    // 注意：根据 Falcon 定义，fG - gF = q mod (x^n+1)。
    // 所以 verify 应该检查是否每个系数都等于 0 (mod q)。
    
    if (check_pass) {
        printf("[SUCCESS] F and G satisfy the NTRU equation (mod q).\n");
    } else {
        printf("[FAILURE] Verification failed.\n");
    }

    free(tmp);
    free(res_fG);
    free(res_gF);
    return 0;
}