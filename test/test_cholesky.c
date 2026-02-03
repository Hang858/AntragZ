#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "common.h"
#include "fft.h"

// 引入 cholesky 函数声明 (或者在头文件中声明)
void run_cholesky(poly *a11, poly *a21, poly *a22,
                  const poly *p11, const poly *p21, const poly *p22);

// ==========================================
// 辅助函数
// ==========================================

// 生成随机多项式
void gen_random_poly(poly *p) {
    for (int i = 0; i < ANTRAG_D; i++) {
        p->coeffs[i] = ((double)(rand() % 2000 - 1000)) / 100.0;
    }
}

// 生成自伴多项式 (f = f*)
// 在时域中满足 f[0]是实数, f[i] = -f[N-i] (针对 X^N+1 环的一种定义，或者 f[i]=conj(f[-i]))
// 注意：对于 X^N+1 环，f* 的定义通常是 f*(x) = f(x^-1) = f(-x^(N-1))
// 具体构造取决于你的 adj 实现。
// 这里我们用一种通用的方法构造正定矩阵：P = B * B*
void gen_positive_definite_matrix(poly *p11, poly *p21, poly *p22) {
    poly b11, b21, b22;
    gen_random_poly(&b11);
    gen_random_poly(&b21);
    gen_random_poly(&b22);

    // 计算 B*
    poly b11_adj, b21_adj, b22_adj;
    b11_adj = b11; poly_adj_fft(&b11_adj, ANTRAG_LOGD);
    b21_adj = b21; poly_adj_fft(&b21_adj, ANTRAG_LOGD);
    b22_adj = b22; poly_adj_fft(&b22_adj, ANTRAG_LOGD);

    // 构造 P = B * B*
    // B = [b11, 0; b21, b22]
    // P11 = b11 * b11*
    // P21 = b21 * b11*
    // P22 = b21 * b21* + b22 * b22*
    
    // 这里需要用到时域或频域乘法，为了方便我们直接用 FFT 乘法接口
    // 注意：输入需要在频域还是时域取决于你的 poly_mul_fft 接口
    // 假设 poly_mul_fft 需要频域输入
    
    FFT(&b11, ANTRAG_LOGD); FFT(&b11_adj, ANTRAG_LOGD);
    FFT(&b21, ANTRAG_LOGD); FFT(&b21_adj, ANTRAG_LOGD);
    FFT(&b22, ANTRAG_LOGD); FFT(&b22_adj, ANTRAG_LOGD);

    poly_mul_fft(p11, &b11, ANTRAG_LOGD);      // p11暂存 b11
    poly_mul_fft(p11, &b11_adj, ANTRAG_LOGD);  // p11 = b11 * b11* (此时在频域)
    
    poly_mul_fft(p21, &b21, ANTRAG_LOGD);      // p21暂存 b21
    poly_mul_fft(p21, &b11_adj, ANTRAG_LOGD);  // p21 = b21 * b11*

    poly t1, t2;
    t1 = b21; poly_mul_fft(&t1, &b21_adj, ANTRAG_LOGD); // b21*b21*
    t2 = b22; poly_mul_fft(&t2, &b22_adj, ANTRAG_LOGD); // b22*b22*
    poly_add(p22, &t1, ANTRAG_LOGD);
    poly_add(p22, &t2, ANTRAG_LOGD); // p22 = sum

    // 转回时域作为 Cholesky 的输入
    invFFT(p11, ANTRAG_LOGD);
    invFFT(p21, ANTRAG_LOGD);
    invFFT(p22, ANTRAG_LOGD);
}

// 计算两个多项式的最大误差
double get_max_diff(const poly *a, const poly *b) {
    double max_diff = 0.0;
    for (int i = 0; i < ANTRAG_D; i++) {
        double diff = fabs(a->coeffs[i] - b->coeffs[i]);
        if (diff > max_diff) max_diff = diff;
    }
    return max_diff;
}

// ==========================================
// 核心测试
// ==========================================
int main() {
    srand(time(NULL));
    printf("=== Cholesky Decomposition Test ===\n");

    poly p11, p21, p22;
    poly a11, a21, a22;
    
    // 1. 生成测试数据 (正定矩阵)
    printf("[1] Generating positive definite matrix P...\n");
    gen_positive_definite_matrix(&p11, &p21, &p22);

    // 2. 运行 Cholesky
    printf("[2] Running Cholesky decomposition...\n");
    run_cholesky(&a11, &a21, &a22, &p11, &p21, &p22);

    // 3. 验证结果: Reconstruct P' = A * A*
    printf("[3] Verifying result (A * A* == P)...\n");
    
    // 转换到频域进行乘法验证
    poly ta11 = a11, ta21 = a21, ta22 = a22;
    FFT(&ta11, ANTRAG_LOGD);
    FFT(&ta21, ANTRAG_LOGD);
    FFT(&ta22, ANTRAG_LOGD);

    poly ta11_adj = a11, ta21_adj = a21, ta22_adj = a22;
    poly_adj_fft(&ta11_adj, ANTRAG_LOGD); FFT(&ta11_adj, ANTRAG_LOGD);
    poly_adj_fft(&ta21_adj, ANTRAG_LOGD); FFT(&ta21_adj, ANTRAG_LOGD);
    poly_adj_fft(&ta22_adj, ANTRAG_LOGD); FFT(&ta22_adj, ANTRAG_LOGD);

    // Calc P'_11 = a11 * a11*
    poly p11_recalc = ta11;
    poly_mul_fft(&p11_recalc, &ta11_adj, ANTRAG_LOGD);
    invFFT(&p11_recalc, ANTRAG_LOGD);

    // Calc P'_21 = a21 * a11*
    poly p21_recalc = ta21;
    poly_mul_fft(&p21_recalc, &ta11_adj, ANTRAG_LOGD);
    invFFT(&p21_recalc, ANTRAG_LOGD);

    // Calc P'_22 = a21*a21* + a22*a22*
    poly t1 = ta21; poly_mul_fft(&t1, &ta21_adj, ANTRAG_LOGD);
    poly t2 = ta22; poly_mul_fft(&t2, &ta22_adj, ANTRAG_LOGD);
    poly p22_recalc;
    // Assuming poly_add works in place or similar
    // We need to clear p22_recalc first or verify poly_add interface
    // Here assume poly_add(dest, src, logn) adds src to dest
    memset(&p22_recalc, 0, sizeof(poly)); // Zero out
    // Re-FFT because poly_add usually works in FFT domain for efficiency or standard add
    // Wait, poly_add is usually element-wise add, works in both domains if size matches.
    // Let's assume it adds coefficients.
    invFFT(&t1, ANTRAG_LOGD);
    invFFT(&t2, ANTRAG_LOGD);
    for(int i=0; i<ANTRAG_D; i++) p22_recalc.coeffs[i] = t1.coeffs[i] + t2.coeffs[i];

    // 4. Check Errors
    double err1 = get_max_diff(&p11, &p11_recalc);
    double err2 = get_max_diff(&p21, &p21_recalc);
    double err3 = get_max_diff(&p22, &p22_recalc);

    printf("  Max Diff P11: %.2e\n", err1);
    printf("  Max Diff P21: %.2e\n", err2);
    printf("  Max Diff P22: %.2e\n", err3);

    if (err1 < 1e-9 && err2 < 1e-9 && err3 < 1e-9) {
        printf("[PASS] Verification successful!\n");
        return 0;
    } else {
        printf("[FAIL] Verification failed (Errors too large).\n");
        return 1;
    }
}