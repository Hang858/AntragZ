
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "common.h"
#include "fft.h"
#include "rng.h"    


static inline double fpr_sqr_scalar(double v) {
    return v * v;
}


static void decode_odd(int8_t u[ANTRAG_D], const poly *utilde) {
    uint8_t umod2 = 0;
    int8_t ui, wi = 0;
    size_t worst_coeff = 0;
    double maxdiff = -1, uitilde, diff;

    for(size_t i=0; i<ANTRAG_D; i++) {
        uitilde = utilde->coeffs[i];
        ui = lrint(uitilde); 
        umod2 ^= ui;             
        diff = fabs(uitilde - (double)ui);
        
        if(diff > maxdiff) {
            worst_coeff = i;
            maxdiff = diff;
            if(uitilde > (double)ui) wi = ui + 1;
            else wi = ui - 1;
        }
        u[i] = ui;
    }
    
    if((umod2 & 1) == 0)
        u[worst_coeff] = wi;
}


int keygen_fg_impl(secret_key_fg *sk) {

    double z[ANTRAG_D/2];
    double af[ANTRAG_D/2], ag[ANTRAG_D/2];
    
    double f_real[ANTRAG_D], g_real[ANTRAG_D];
    

    const double alow  = 0.5*(ANTRAG_ALPHA + 1.0/ANTRAG_ALPHA) - 0.5*ANTRAG_XI*(ANTRAG_ALPHA-1.0/ANTRAG_ALPHA);
    const double ahigh = 0.5*(ANTRAG_ALPHA + 1.0/ANTRAG_ALPHA) + 0.5*ANTRAG_XI*(ANTRAG_ALPHA-1.0/ANTRAG_ALPHA);
    const double qlow  = ((double)ANTRAG_Q)*alow*alow;    // r^2
    const double qhigh = ((double)ANTRAG_Q)*ahigh*ahigh;  // R^2
    
    const double qlow2 = ((double)ANTRAG_Q)/(ANTRAG_ALPHA*ANTRAG_ALPHA);  // q / alpha^2
    const double qhigh2= ((double)ANTRAG_Q)*ANTRAG_ALPHA*ANTRAG_ALPHA;    // q * alpha^2

    double *r = (double*)malloc(2 * ANTRAG_D * sizeof(double));
    if (!r) return -1;

    int trials = 0;
    int check;

    do {
        trials++;
        
        sample_fpr(r, 2 * ANTRAG_D);

        for(size_t i=0; i<ANTRAG_D/2; i++) {

            z[i] = sqrt(qlow + (qhigh - qlow)*r[i]); // row
            
            double theta1 = M_PI/2 * r[i + ANTRAG_D/2];
            af[i] = z[i] * cos(theta1);
            ag[i] = z[i] * sin(theta1);  // (x, y) -> (r cosθ, r sinθ)
            
            double theta2 = 2*M_PI * r[i + 2*ANTRAG_D/2];
            f_real[i] = af[i] * cos(theta2); // 实部
            f_real[i+ANTRAG_D/2] = af[i] * sin(theta2); // 虚部  //x * e^(i theta2)
            
            double theta3 = 2*M_PI * r[i + 3*ANTRAG_D/2];
            g_real[i] = ag[i] * cos(theta3); // 实部
            g_real[i+ANTRAG_D/2] = ag[i] * sin(theta3); // 虚部  //y * e^(i theta3)
        }

        for(size_t i=0; i<ANTRAG_D; i++) {
            sk->b10.coeffs[i] = f_real[i];
            sk->b11.coeffs[i] = g_real[i];
        }

        invFFT(&sk->b10, ANTRAG_LOGD);
        invFFT(&sk->b11, ANTRAG_LOGD);


        decode_odd(sk->f, &(sk->b10));
        decode_odd(sk->g, &(sk->b11));

        for(size_t i=0; i<ANTRAG_D; i++) {
            sk->b10.coeffs[i] = (double)sk->f[i];
            sk->b11.coeffs[i] = (double)sk->g[i];
        }

        FFT(&sk->b10, ANTRAG_LOGD);
        FFT(&sk->b11, ANTRAG_LOGD);

        check = 0;
        for(size_t i=0; i<ANTRAG_D/2; i++) {

            double norm_f = fpr_sqr_scalar(sk->b10.coeffs[i]) + 
                            fpr_sqr_scalar(sk->b10.coeffs[i+ANTRAG_D/2]);
            
            double norm_g = fpr_sqr_scalar(sk->b11.coeffs[i]) + 
                            fpr_sqr_scalar(sk->b11.coeffs[i+ANTRAG_D/2]);
            
            double zi = norm_f + norm_g;
            
            if(zi < qlow2 || zi > qhigh2) {
                check = 1;
                break;
            }
        }

    } while(check);

    free(r);
    return trials;
}