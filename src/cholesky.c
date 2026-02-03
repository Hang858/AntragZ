#include <complex.h>
#include "common.h"
#include "fft.h"

// ---------------------------------------------------------
// 辅助函数：手动处理 Split Layout (Re: 0..N/2-1, Im: N/2..N-1)
// ---------------------------------------------------------

static void poly_div_fft(poly *out, const poly *a, const poly *b) {
    int n = ANTRAG_D / 2; // 复数对的个数 (256)
    
    for (int i = 0; i < n; i++) {
        // 1. 从 Split Layout 构造复数
        double a_re = a->coeffs[i];
        double a_im = a->coeffs[i + n];
        double complex ca = a_re + I * a_im;

        double b_re = b->coeffs[i];
        double b_im = b->coeffs[i + n];
        double complex cb = b_re + I * b_im;

        // 2. 执行复数除法
        // 注意：cb 接近 0 时可能会产生巨大数值，但 P 正定性应保证 cb 模长足够大
        double complex res = ca / cb;

        // 3. 存回 Split Layout
        out->coeffs[i]     = creal(res);
        out->coeffs[i + n] = cimag(res);
    }
}

static void poly_sqrt_fft(poly *out, const poly *src) {
    int n = ANTRAG_D / 2;
    
    for (int i = 0; i < n; i++) {
        // 1. 从 Split Layout 构造复数
        double s_re = src->coeffs[i];
        double s_im = src->coeffs[i + n];
        double complex cs = s_re + I * s_im;

        // 2. 执行复数开方
        double complex res = csqrt(cs);

        // 3. 存回 Split Layout
        out->coeffs[i]     = creal(res);
        out->coeffs[i + n] = cimag(res);
    }
}


/**
 * @brief 执行 Cholesky 分解 P = A * A^T
 * * @param a11 [输出] 结果矩阵 A 的 (1,1) 元素
 * @param a21 [输出] 结果矩阵 A 的 (2,1) 元素
 * @param a22 [输出] 结果矩阵 A 的 (2,2) 元素
 * @param p11 [输入] 输入矩阵 P 的 (1,1) 元素 (自伴)
 * @param p21 [输入] 输入矩阵 P 的 (2,1) 元素
 * @param p22 [输入] 输入矩阵 P 的 (2,2) 元素 (自伴)
 */
void run_cholesky(poly *a11, poly *a21, poly *a22,
                  const poly *p11, const poly *p21, const poly *p22) 
{

    poly t_p11 = *p11;
    poly t_p21 = *p21;
    poly t_p22 = *p22;
    
    FFT(&t_p11, ANTRAG_LOGD);
    FFT(&t_p21, ANTRAG_LOGD);
    FFT(&t_p22, ANTRAG_LOGD);

    poly_sqrt_fft(a11, &t_p11);

    poly_div_fft(a21, &t_p21, a11);

    poly t_schur_term;

    poly_mul_autoadj_fft(&t_schur_term, a21, ANTRAG_LOGD);

    poly_sub(&t_p22, &t_schur_term, ANTRAG_LOGD);

    poly_sqrt_fft(a22, &t_p22);

    invFFT(a11, ANTRAG_LOGD);
    invFFT(a21, ANTRAG_LOGD);
    invFFT(a22, ANTRAG_LOGD);
}