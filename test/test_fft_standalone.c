#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h> // 需要安装 FFTW 开发库
#include "common.h"
#include "fft.h"   // 这里声明的是你的 KISS FFT 版本接口

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// =========================================================
// 参考实现 (Reference): 原版 FFTW 的 invFFT 代码
// =========================================================
void invFFT_FFTW(poly *p, int logn) {
   int n = 1 << logn;
    fftw_complex *in, *out;
    fftw_plan plan;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    
    for (int k = 0; k < n / 2; k++) {
        double re = p->coeffs[k];
        double im = p->coeffs[k + n/2];

        in[k][0] = re;
        in[k][1] = im;
        
        int sym_k = n - 1 - k;
        in[sym_k][0] = re;
        in[sym_k][1] = -im;
    }

    plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // 3. Untwist & Normalize
    // b_j = out[j] * psi^(-j) / n
    for (int j = 0; j < n; j++) {
        double re = out[j][0];
        double im = out[j][1];
        
        double angle = -1.0 * (M_PI * j) / n;
        double cos_a = cos(angle);
        double sin_a = sin(angle);
        
        double val_re = re * cos_a - im * sin_a;
        
        p->coeffs[j] = val_re / n;
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}

// =========================================================
// 主测试逻辑
// =========================================================
int main() {
    printf("=== Inverse FFT (invFFT) Comparison (KISS vs FFTW) ===\n");

    int logn = 9;
    int n = 1 << logn; // 512
    poly p_kiss, p_fftw;

    // 1. 构造测试输入
    // 输入是频域数据（前半部分实部，后半部分虚部）
    for (int i = 0; i < n; i++) {
        // 使用一些非对称的数据来测试复数运算
        if (i < n/2) {
            p_kiss.coeffs[i] = (double)i * 0.5; // 实部
        } else {
            p_kiss.coeffs[i] = (double)(i - n/2) * -0.5; // 虚部
        }
        
        // 制造一点脉冲
        if (i == 10) p_kiss.coeffs[i] += 100.0;
        
        p_fftw.coeffs[i] = p_kiss.coeffs[i];
    }

    // 2. 分别执行逆变换
    // KISS 实现
    invFFT(&p_kiss, logn);
    
    // FFTW 参考实现
    invFFT_FFTW(&p_fftw, logn);

    // 3. 逐项对比时域系数
    printf("\nComparing Time Domain Coefficients:\n");
    printf("Index | KISS Value             | FFTW Value             | Diff\n");
    printf("------+------------------------+------------------------+--------\n");

    int diff_count = 0;
    double max_diff = 0.0;

    for (int k = 0; k < n; k++) {
        double val_kiss = p_kiss.coeffs[k];
        double val_fftw = p_fftw.coeffs[k];

        double diff = fabs(val_kiss - val_fftw);
        
        if (diff > 1e-8) {
            if (diff_count < 10) { 
                printf("%3d   | %22.10f | %22.10f | %.1e <--- FAIL\n",
                       k, val_kiss, val_fftw, diff);
            }
            if (diff > max_diff) max_diff = diff;
            diff_count++;
        }
    }

    printf("------+------------------------+------------------------+--------\n");
    if (diff_count == 0) {
        printf("[PASS] Inverse FFT implementations match perfectly!\n");
        return 0;
    } else {
        printf("[FAIL] Found %d mismatches. Max diff: %e\n", diff_count, max_diff);
        return 1;
    }
}