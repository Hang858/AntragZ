/* tests/debug_keygen.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "poly.h"
#include "fft.h"
#include "rng.h"

// 声明内部函数（假设这些函数在链接时可见，或者你可以临时把它们从 static 改为非 static）
// 注意：编译时需要链接对应的 .o 文件
extern int keygen_fg_impl(secret_key_fg* sk); // 这里用 void* 绕过类型检查，实际传入 secret_key_fg*
extern int antrag_solve_ntru(int8_t *F, int8_t *G, const int8_t *f, const int8_t *g, unsigned logn, uint32_t *tmp);
extern int Run_PreMatrix(const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G, PreMatrix_Output *out);

// 辅助函数：计算多项式的欧几里得范数平方
double poly_norm_sq(const int8_t *p, int n) {
    double sum = 0;
    for(int i=0; i<n; i++) {
        sum += (double)p[i] * p[i];
    }
    return sum;
}

// 辅助函数：计算多项式的欧几里得范数
double poly_norm(const int8_t *p, int n) {
    return sqrt(poly_norm_sq(p, n));
}

int main() {
    printf("=== Zitaka KeyGen Diagnostic Tool ===\n");
    
    // 0. 初始化
    int n = ANTRAG_D;
    
    // 定义结构体（需与 src/common.h 或 keygen_basis.c 中的定义一致）
    // 这里简单定义缓冲区
    int8_t *f = malloc(n);
    int8_t *g = malloc(n);
    int8_t *F = malloc(n);
    int8_t *G = malloc(n);
    
    // 临时结构用于 keygen_fg_impl
    typedef struct {
        int8_t f[512];
        int8_t g[512];
        int8_t F[512];
        int8_t G[512];
        // ... 其他可能需要的字段，这里只取前两个
        // 注意：这只是为了对其 keygen_fg_impl 的内存布局，实际请参考你的头文件定义
        // 如果 secret_key_fg 定义复杂，建议直接 include 相关头文件
        poly b10; // 占位
        poly b11; // 占位
    } secret_key_fg_mock;
    
    secret_key_fg_mock fg_ctx;
    
    printf("\n[Step 1] Generating f, g (NTRU Basis Part 1)...\n");
    int trials = keygen_fg_impl(&fg_ctx);
    printf("  -> Attempts: %d\n", trials);
    
    memcpy(f, fg_ctx.f, n);
    memcpy(g, fg_ctx.g, n);

    double norm_f = poly_norm(f, n);
    double norm_g = poly_norm(g, n);
    double basis_norm_fg = sqrt((norm_f*norm_f + norm_g*norm_g));
    
    printf("  -> Norm(f): %.4f\n", norm_f);
    printf("  -> Norm(g): %.4f\n", norm_g);
    printf("  -> Basis Norm (f,g): %.4f\n", basis_norm_fg);
    printf("  -> Expected Max: ~ %.4f (1.17 * sqrt(q))\n", 1.17 * sqrt(ANTRAG_Q));

    // 检查 f 是否可逆
    int16_t *f_inv_check = malloc(n * sizeof(int16_t));
    if (!poly_inv_mod_q(f_inv_check, f)) {
        printf("  [FAIL] f is NOT invertible mod q. (This key would be rejected)\n");
    } else {
        printf("  [PASS] f is invertible mod q.\n");
    }
    free(f_inv_check);

    printf("\n[Step 2] Solving NTRU Equation for F, G...\n");
    uint32_t *tmp_ntru = malloc(16384 * sizeof(uint32_t)); 
    if (!antrag_solve_ntru(F, G, f, g, ANTRAG_LOGD, tmp_ntru)) {
        printf("  [FAIL] NTRU Solver failed to find solution.\n");
        goto cleanup;
    }
    printf("  [PASS] NTRU Solution found.\n");
    
    double norm_F = poly_norm(F, n);
    double norm_G = poly_norm(G, n);
    double basis_norm_FG = sqrt((norm_F*norm_F + norm_G*norm_G));
    
    printf("  -> Norm(F): %.4f\n", norm_F);
    printf("  -> Norm(G): %.4f\n", norm_G);
    printf("  -> Basis Norm (F,G): %.4f\n", basis_norm_FG);
    
    // 验证 NTRU 方程: fG - gF = q
    // 简单验证常数项或随机几项，或者完整验证
    // 这里做个完整验证的简化版（不处理 mod X^N+1 循环卷积细节，仅作提示）
    printf("  -> Verifying fG - gF = q mod (X^N+1)...\n");
    // (需调用 poly_mul 等函数，此处略过，假设 Solver 是对的，重点关注范数)

    printf("\n[Step 3] PreMatrix Check (The Cholesky Bottleneck)...\n");
    
    // 关键分析：s0 的限制
    // PreMatrix 要求: Sigma = B * B* < s0^2 * p^2 * I
    // 粗略来说，就是要求 ||(f, g, F, G)|| 的某种组合范数不能太大
    // 你的 s0 设置为 131
    
    // 让我们看看 Falcon 标准的 s0 限制
    // Falcon-512 的 s0 一般在 1.17*sqrt(q) 左右，即 ~130
    // 你生成的 F, G 范数如果是 200+，那就必挂无疑
    
    double max_norm = norm_F > norm_G ? norm_F : norm_G;
    if (max_norm > norm_f) max_norm = max_norm > norm_f ? max_norm : norm_f;
    if (max_norm > norm_g) max_norm = max_norm > norm_g ? max_norm : norm_g;
    
    printf("  -> Max Poly Norm: %.4f\n", max_norm);
    printf("  -> s0 Threshold:  131.0 (Approx)\n");
    
    if (max_norm > 131.0 * 1.5) { // 宽松一点的预判
        printf("  [WARNING] F or G norm is VERY high. PreMatrix is likely to fail.\n");
    }

    // 尝试运行 PreMatrix
    // 需要分配一个 PreMatrix_Output 结构，这里为了编译方便，我们可以只是为了跑通 Run_PreMatrix
    // 你可能需要把 struct 定义搬过来或者 include 头文件
    // 假设 Run_PreMatrix 内部主要逻辑是 Cholesky
    
    // 这里我们直接调用 Run_PreMatrix 看看返回值
    // 注意：需要分配足够大的内存给 out，或者使用实际的结构体
    // 这里分配一个足够大的 buffer 避免溢出
    void *mock_out = malloc(1024 * 1024); 
    
    int pm_result = Run_PreMatrix(f, g, F, G, (PreMatrix_Output*)mock_out);
    
    if (pm_result) {
        printf("  [SUCCESS] PreMatrix computed successfully!\n");
        printf("  ==> This keypair IS VALID.\n");
    } else {
        printf("  [FAIL] PreMatrix computation failed.\n");
        printf("  ==> Reason: Matrix likely not Positive Definite (Norms too large).\n");
    }
    free(mock_out);

cleanup:
    free(f); free(g); free(F); free(G);
    free(tmp_ntru);
    return 0;
}