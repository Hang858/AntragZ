#include <complex.h>
#include "common.h"
#include "fft.h"



// static void poly_div_fft(poly *out, const poly *a, const poly *b) {
//     int n = ANTRAG_D / 2; // 复数对的个数 (256)
    
//     for (int i = 0; i < n; i++) {
//         // 1. 从 Split Layout 构造复数
//         double a_re = a->coeffs[i];
//         double a_im = a->coeffs[i + n];
//         double complex ca = a_re + I * a_im;

//         double b_re = b->coeffs[i];
//         double b_im = b->coeffs[i + n];
//         double complex cb = b_re + I * b_im;

//         // 2. 执行复数除法
//         // 注意：cb 接近 0 时可能会产生巨大数值，但 P 正定性应保证 cb 模长足够大
//         double complex res = ca / cb;

//         // 3. 存回 Split Layout
//         out->coeffs[i]     = creal(res);
//         out->coeffs[i + n] = cimag(res);
//     }
// }

// static void poly_sqrt_fft(poly *out, const poly *src) {
//     int n = ANTRAG_D / 2;
    
//     for (int i = 0; i < n; i++) {
//         // 1. 从 Split Layout 构造复数
//         double s_re = src->coeffs[i];
//         double s_im = src->coeffs[i + n];
//         double complex cs = s_re + I * s_im;

//         // 2. 执行复数开方
//         double complex res = csqrt(cs);

//         // 3. 存回 Split Layout
//         out->coeffs[i]     = creal(res);
//         out->coeffs[i + n] = cimag(res);
//     }
// }



#include <math.h>
#include <stdio.h>
#include "fft.h"

// 设定一个极小的安全阈值，防止因为 0.00000000000001 变成负数导致 NaN
#define EPSILON 1e-15

void run_cholesky(poly *a11, poly *a21, poly *a22,
                 const poly *p11, const poly *p21, const poly *p22) 
{
    int n = ANTRAG_D;
    int hn = n / 2;

    // 临时拷贝输入并转到 FFT 域
    poly tp11 = *p11, tp21 = *p21, tp22 = *p22;
    FFT(&tp11, ANTRAG_LOGD);
    FFT(&tp21, ANTRAG_LOGD);
    FFT(&tp22, ANTRAG_LOGD);

    for (int i = 0; i < hn; i++) {
        // --- 1. 计算 a11 = sqrt(p11) ---
        // p11 是正定矩阵对角线，点值必为正实数
        double p11_re = tp11.coeffs[i];
        if (p11_re < EPSILON) p11_re = EPSILON; 
        double val_a11 = sqrt(p11_re);
        
        a11->coeffs[i] = val_a11;
        a11->coeffs[i + hn] = 0.0; // 对角线元素 FFT 后虚部恒为 0

        // --- 2. 计算 a21 = p21 / a11 ---
        // 注意：由于 a11 是实数点值，直接除即可
        double p21_re = tp21.coeffs[i];
        double p21_im = tp21.coeffs[i + hn];
        
        a21->coeffs[i] = p21_re / val_a11;
        a21->coeffs[i + hn] = p21_im / val_a11;

        // --- 3. 计算 a22 = sqrt(p22 - |a21|^2) ---
        double a21_re = a21->coeffs[i];
        double a21_im = a21->coeffs[i + hn];
        double a21_mag_sq = a21_re * a21_re + a21_im * a21_im;
        
        double p22_schur = tp22.coeffs[i] - a21_mag_sq;
        if (p22_schur < EPSILON) p22_schur = EPSILON; // 再次防止负值产生的 NaN
        
        a22->coeffs[i] = sqrt(p22_schur);
        a22->coeffs[i + hn] = 0.0; // 对角线虚部恒为 0
    }

    // 转回系数域
    invFFT(a11, ANTRAG_LOGD);
    invFFT(a21, ANTRAG_LOGD);
    invFFT(a22, ANTRAG_LOGD);
}