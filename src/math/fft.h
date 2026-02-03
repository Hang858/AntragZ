#ifndef FFT_H
#define FFT_H

#include "common.h"

void FFT(poly *p, int logn);

void invFFT(poly *p, int logn);

void poly_mul_fft(poly *a, const poly *b, int logn);

void poly_add(poly *a, const poly *b, int logn);

void poly_sub(poly *a, const poly *b, int logn);

void poly_adj_fft(poly *a, int logn);
void poly_muladj_fft(poly *a, const poly *b, int logn);

void poly_mulselfadj_fft(poly *a, int logn);

void poly_invnorm2_fft(poly *d, const poly *a, const poly *b, int logn);

void poly_mul_autoadj_fft(poly *a, const poly *b, unsigned logn);

void poly_add_muladj_fft(poly *d, const poly *F, const poly *G, 
                         const poly *f, const poly *g, unsigned logn);

#endif