#include <fftw3.h>
#include <math.h>
#include "common.h"
#define ANTRAG_D 512

void FFT(poly *p, int logn) {
    int n = 1 << logn;
    fftw_complex *in, *out;
    fftw_plan plan;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);


    for (int j = 0; j < n; j++) {
        double re = p->coeffs[j]; // 输入是实数
        double angle = (M_PI * j) / n;
        
        // twist = cos + i*sin
        // in[j] = re * twist
        in[j][0] = re * cos(angle);
        in[j][1] = re * sin(angle);
    }

    plan = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    for (int k = 0; k < n / 2; k++) {
        p->coeffs[k] = out[k][0];           // 实部
        p->coeffs[k + n/2] = out[k][1];     // 虚部
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}


void invFFT(poly *p, int logn) {
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



void poly_mul_fft(poly *a, const poly *b, int logn) {
    int n = 1 << logn;
    int half = n / 2;
    for (int i = 0; i < half; i++) {

        double a_re = a->coeffs[i];
        double a_im = a->coeffs[i + half];
        double b_re = b->coeffs[i];
        double b_im = b->coeffs[i + half];

        a->coeffs[i] = a_re * b_re - a_im * b_im;
        a->coeffs[i + half] = a_re * b_im + a_im * b_re;
    }
}

static void fpr_mul_complex(double *re, double *im, double a_re, double a_im, double b_re, double b_im) {
    *re = a_re * b_re - a_im * b_im;
    *im = a_re * b_im + a_im * b_re;
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
        a->coeffs[i] = -a->coeffs[i]; // 虚部取反
    }
}

void poly_muladj_fft(poly *a, const poly *b, int logn) {
    int n = 1 << logn;
    int hn = n / 2;
    for (int i = 0; i < hn; i++) {
        double a_re = a->coeffs[i];
        double a_im = a->coeffs[i + hn];
        double b_re = b->coeffs[i];
        double b_im = -b->coeffs[i + hn]; // b 的虚部取反
        
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


