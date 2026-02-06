#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "fft.h"

// 1. 配置 KISS FFT
#define kiss_fft_scalar double 
#include "kiss_fft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// =========================================================
// 2. 动态配置管理 (Fix for dynamic logn)
// =========================================================

// 最大支持 logn = 10 (N=1024)，足够 Zitaka/Falcon 使用
#define MAX_LOGN 11 

// 使用数组缓存每个 logn 对应的配置
// cache[9] 存 logn=9 的配置，cache[8] 存 logn=8 的配置...
static kiss_fft_cfg cfg_fwd_cache[MAX_LOGN] = {0}; 
static kiss_fft_cfg cfg_inv_cache[MAX_LOGN] = {0};

// 静态缓冲区 (最大 N=512)
static kiss_fft_cpx fft_in[ANTRAG_D];
static kiss_fft_cpx fft_out[ANTRAG_D];

static void ensure_fft_initialized(int logn) {
    if (logn < 0 || logn >= MAX_LOGN) {
        // 异常保护：超出范围
        return; 
    }

    int n = 1 << logn; 
    
    // 检查当前 logn 对应的配置是否已分配
    // 如果这一层 logn 还没初始化过，就分配它
    if (!cfg_fwd_cache[logn]) {
        cfg_fwd_cache[logn] = kiss_fft_alloc(n, 0, NULL, NULL);
    }
    if (!cfg_inv_cache[logn]) {
        cfg_inv_cache[logn] = kiss_fft_alloc(n, 1, NULL, NULL);
    }
}

// =========================================================
// 3. 核心 FFT/invFFT
// =========================================================

void FFT(poly *p, int logn) {
    // 1. 确保当前维度的配置已就绪
    ensure_fft_initialized(logn);
    
    // 获取当前维度的配置
    kiss_fft_cfg cfg = cfg_inv_cache[logn]; // 注意：用 inv 配置 (指数+1)

    int n = 1 << logn;
    
    // 预处理 (Twist)
    for (int j = 0; j < n; j++) {
        double re = p->coeffs[j]; 
        double angle = (M_PI * j) / n;
        
        double twist_re = cos(angle);
        double twist_im = sin(angle);

        fft_in[j].r = re * twist_re;
        fft_in[j].i = re * twist_im;
    }

    // 执行变换
    kiss_fft(cfg, fft_in, fft_out);

    // 后处理 (Split)
    int half_n = n >> 1;
    for (int k = 0; k < half_n; k++) {
        p->coeffs[k]          = fft_out[k].r; 
        p->coeffs[k + half_n] = fft_out[k].i; 
    }
}


void invFFT(poly *p, int logn) {
    ensure_fft_initialized(logn);
    
    // 获取当前维度的配置
    kiss_fft_cfg cfg = cfg_fwd_cache[logn]; // 注意：用 fwd 配置 (指数-1)

    int n = 1 << logn;
    int half_n = n >> 1;

    // 组装输入
    for (int k = 0; k < half_n; k++) {
        double re = p->coeffs[k];
        double im = p->coeffs[k + half_n];

        fft_in[k].r = re;
        fft_in[k].i = im;
        
        int sym_k = n - 1 - k;
        fft_in[sym_k].r = re;
        fft_in[sym_k].i = -im;
    }

    // 执行变换
    kiss_fft(cfg, fft_in, fft_out);

    // 后处理 (Untwist & Normalize)
    for (int j = 0; j < n; j++) {
        double re = fft_out[j].r;
        double im = fft_out[j].i;
        
        double angle = -1.0 * (M_PI * j) / n;
        double cos_a = cos(angle);
        double sin_a = sin(angle);
        
        double val_re = re * cos_a - im * sin_a;
        p->coeffs[j] = val_re / n;
    }
}


static void fpr_mul_complex(double *re, double *im, double a_re, double a_im, double b_re, double b_im) {
    *re = a_re * b_re - a_im * b_im;
    *im = a_re * b_im + a_im * b_re;
}

void poly_mul_fft(poly *a, const poly *b, int logn) {
    int n = 1 << logn;
    int half = n / 2;
    for (int i = 0; i < half; i++) {
        double a_re = a->coeffs[i];
        double a_im = a->coeffs[i + half];
        double b_re = b->coeffs[i];
        double b_im = b->coeffs[i + half];

        a->coeffs[i]        = a_re * b_re - a_im * b_im;
        a->coeffs[i + half] = a_re * b_im + a_im * b_re;
    }
}

void poly_add(poly *a, const poly *b, int logn) {
    int n = 1 << logn;
    for (int i = 0; i < n; i++) {
        a->coeffs[i] += b->coeffs[i];
    }
}

void poly_sub(poly *a, const poly *b, int logn) {
    int n = 1 << logn;
    for (int i = 0; i < n; i++) {
        a->coeffs[i] -= b->coeffs[i];
    }
}

void poly_adj_fft(poly *a, int logn) {
    int n = 1 << logn;
    int hn = n / 2;
    for (int i = hn; i < n; i++) {
        a->coeffs[i] = -a->coeffs[i]; 
    }
}

void poly_muladj_fft(poly *a, const poly *b, int logn) {
    int n = 1 << logn;
    int hn = n / 2;
    for (int i = 0; i < hn; i++) {
        double a_re = a->coeffs[i];
        double a_im = a->coeffs[i + hn];
        double b_re = b->coeffs[i];
        double b_im = -b->coeffs[i + hn]; 
        
        double res_re, res_im;
        fpr_mul_complex(&res_re, &res_im, a_re, a_im, b_re, b_im);
        
        a->coeffs[i] = res_re;
        a->coeffs[i + hn] = res_im;
    }
}

void poly_mulselfadj_fft(poly *a, int logn) {
    int n = 1 << logn;
    int hn = n / 2;
    for (int i = 0; i < hn; i++) {
        double re = a->coeffs[i];
        double im = a->coeffs[i + hn];
        a->coeffs[i] = re * re + im * im;
        a->coeffs[i + hn] = 0.0;
    }
}

void poly_invnorm2_fft(poly *d, const poly *a, const poly *b, int logn) {
    int n = 1 << logn;
    int hn = n / 2;
    for (int i = 0; i < hn; i++) {
        double a_re = a->coeffs[i];
        double a_im = a->coeffs[i + hn];
        double b_re = b->coeffs[i];
        double b_im = b->coeffs[i + hn];
        
        double val = (a_re*a_re + a_im*a_im) + (b_re*b_re + b_im*b_im);
        d->coeffs[i] = 1.0 / val; 
        d->coeffs[i + hn] = 0.0;
    }
}

void poly_mul_autoadj_fft(poly *a, const poly *b, unsigned logn) {
    size_t n = (size_t)1 << logn;
    size_t hn = n >> 1;
    for (size_t i = 0; i < hn; i++) {
        double b_re = b->coeffs[i]; 
        a->coeffs[i]    *= b_re;
        a->coeffs[i+hn] *= b_re;
    }
}

void poly_add_muladj_fft(poly *d, const poly *F, const poly *G, 
                         const poly *f, const poly *g, unsigned logn) {
    size_t n = (size_t)1 << logn;
    size_t hn = n >> 1;
    for (size_t i = 0; i < hn; i++) {
        double F_re = F->coeffs[i],    F_im = F->coeffs[i+hn];
        double G_re = G->coeffs[i],    G_im = G->coeffs[i+hn];
        double f_re = f->coeffs[i],    f_im = f->coeffs[i+hn];
        double g_re = g->coeffs[i],    g_im = g->coeffs[i+hn];

        double re1 = F_re * f_re + F_im * f_im;
        double im1 = F_im * f_re - F_re * f_im;
        
        double re2 = G_re * g_re + G_im * g_im;
        double im2 = G_im * g_re - G_re * g_im;

        d->coeffs[i]    = re1 + re2;
        d->coeffs[i+hn] = im1 + im2;
    }
}