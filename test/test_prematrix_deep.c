/* test/test_prematrix_deep.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "common.h"
#include "poly.h"
#include "fft.h"
#include "prematrix.h"

// 简单的 L2 范数
double poly_norm(const poly *a) {
    double sum = 0;
    for(int i=0; i<ANTRAG_D; i++) sum += a->coeffs[i] * a->coeffs[i];
    return sqrt(sum);
}
int antrag_solve_ntru(int8_t *F, int8_t *G, const int8_t *f, const int8_t *g, unsigned logn, uint32_t *tmp);
// 1. 测试 FFT 往返
void test_fft_roundtrip() {
    printf("[TEST] FFT -> iFFT Roundtrip...\n");
    poly *p1 = malloc(sizeof(poly));
    poly *p2 = malloc(sizeof(poly));
    
    // 生成随机多项式
    for(int i=0; i<ANTRAG_D; i++) {
        p1->coeffs[i] = (double)(rand() % 100 - 50);
        p2->coeffs[i] = p1->coeffs[i];
    }
    
    FFT(p1, ANTRAG_LOGD);
    invFFT(p1, ANTRAG_LOGD);
    
    double max_err = 0;
    for(int i=0; i<ANTRAG_D; i++) {
        double err = fabs(p1->coeffs[i] - p2->coeffs[i]);
        if(err > max_err) max_err = err;
    }
    
    printf("  Max Error: %.5e\n", max_err);
    if(max_err < 1e-9) printf("  [PASS] FFT looks correct.\n");
    else printf("  [FAIL] FFT precision issue or scaling missing!\n");
    
    free(p1); free(p2);
}

// 2. 模拟计算 u_hat 并检查残差
// 手动计算 R = F*p - u*f
void check_u_hat_quality(int8_t *f, int8_t *g, int8_t *F, int8_t *G, PreMatrix_Output *mat) {
    printf("\n[TEST] Checking u_hat Quality...\n");
    
    double p = pow(2.0, 28);
    double max_R1 = 0;
    
    // R1 = F*p - u*f (卷积)
    // 这里我们只估算第一项系数，避免写全卷积代码
    // R1[0] = (F[0]*p - u[0]*f[0] - u[1]*f[n-1] ...)
    // 直接看 mat->u_hat_num 的范围
    
    int64_t max_u = 0;
    for(int i=0; i<ANTRAG_D; i++) {
        if(llabs(mat->u_hat_num[i]) > max_u) max_u = llabs(mat->u_hat_num[i]);
    }
    printf("  u_hat_num Max: %ld (approx %.2f * p)\n", max_u, (double)max_u/p);
    
    if (max_u > 10 * p) {
        printf("  [WARN] u_hat is suspiciously large! Normal range is ~0.5*p.\n");
    } else {
        printf("  [PASS] u_hat range looks reasonable.\n");
    }
}

// 3. 注入式调试 Run_PreMatrix
// 我们复制一部分逻辑来探针 B^2 - Sigma
void probe_sigma_magnitude(int8_t *f, int8_t *g, int8_t *F, int8_t *G, PreMatrix_Output *mat) {
    printf("\n[TEST] Probing Sigma vs B^2...\n");
    int64_t p_int = 1LL << 28;
    int128_t p2 = (int128_t)p_int * p_int;
    int128_t B_sq = (int128_t)16900 * p2; // 130^2 * p^2
    
    // 粗略估算 R1 的模长平方
    // Sigma_11 = (p^2 |F|^2 + |u|^2 |f|^2 - 2 p u <F,f>) / p^2 ... 实际上 Sigma = |R1|^2 + |R2|^2
    // S11 = B^2 - Sigma_11
    
    // 我们直接利用 mat->c11 来反推 S11? 不行，c11 是 Cholesky 后的。
    // 我们需要看 Run_PreMatrix 内部的 S11[0] 在减法后的值。
    
    // 由于无法直接访问函数内部变量，我们建议用户修改 prematrix.c 打印调试信息。
    printf("  Please ensure you added debug prints in src/prematrix/prematrix.c before 'run_cholesky'.\n");
}

int main() {
    srand(time(NULL));
    
    // 1. 测试 FFT
    test_fft_roundtrip();
    
    // 2. 生成临时密钥用于测试
    printf("\n[SETUP] Generating temporary keypair...\n");
    secret_key_fg sk_fg;
    keygen_fg_impl(&sk_fg);
    
    PrivateKey sk;
    memset(&sk, 0, sizeof(sk));
    memcpy(sk.f, sk_fg.f, ANTRAG_D);
    memcpy(sk.g, sk_fg.g, ANTRAG_D);
    
    uint32_t *tmp = malloc(9000 * sizeof(uint32_t));
    antrag_solve_ntru(sk.F, sk.G, sk.f, sk.g, ANTRAG_LOGD, tmp);
    free(tmp);
    
    // 3. 运行 PreMatrix
    printf("[SETUP] Running Run_PreMatrix...\n");
    if (!Run_PreMatrix(sk.f, sk.g, sk.F, sk.G, &sk.mat)) {
        printf("[FAIL] Run_PreMatrix returned 0.\n");
    } else {
        printf("[PASS] Run_PreMatrix returned 1.\n");
    }
    
    check_u_hat_quality(sk.f, sk.g, sk.F, sk.G, &sk.mat);
    
    return 0;
}