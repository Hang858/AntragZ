#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stddef.h>

#define ANTRAG_D 512        // n
#define ANTRAG_LOGD 9    // log2(n)
#define ANTRAG_Q 12289      // q
#define ANTRAG_ALPHA 1.17   // alpha
#define ANTRAG_ALPHAEPS 0.0
#define ANTRAG_XI 1.0 / 3.0      // xi, 通常接近1
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEPTH_INT_FG   4
#define fpr_add(a, b) ((a) + (b))
#define fpr_sub(a, b) ((a) - (b))
#define fpr_mul(a, b) ((a) * (b))
#define fpr_sqr(a)    ((a) * (a))
#define fpr_inv(a)    (1.0 / (a))
#define fpr_lt(a, b)  ((a) < (b))
#define fpr_rint(a)   rint(a)
#define fpr_of(a)     ((double)(a))
#define fpr_zero      0.0
#define fpr_one       1.0
#define fpr_two       2.0
#define fpr_onehalf   0.5
#define fpr_ptwo31    2147483648.0
#define fpr_ptwo31m1  2147483647.0
#define fpr_mtwo31m1  -2147483647.0
#define fpr_ptwo63m1  9223372036854775807.0
#define fpr_mtwo63m1  -9223372036854775807.0



#define N_MAX 512
#define K_VAL 3
#define L_VAL 9
#define OUTPUT_ROWS 28 // 3*(9-1) + 4
#define MIGD_ROWS 28 // 3*(9-1) + 4
#define COMPRESSED_POOL_SIZE 1024 
#define STACK_BUF_SIZE 512
#define WORK_BUF_SIZE 1600

#define MKN(logn) ((size_t)1 << (logn))

typedef struct {
    double coeffs[ANTRAG_D];
} poly;   //4KB

typedef struct {
    int64_t *coeffs;
    size_t n; // 维度 (2^k)
} poly_int;

typedef struct {
    int8_t f[ANTRAG_D];  // 最终的整数 f
    int8_t g[ANTRAG_D];  // 最终的整数 g
    poly b10;            // 临时变量 / 频域 f
    poly b11;            // 临时变量 / 频域 g
} secret_key_fg;

typedef struct {
    uint32_t p;
    uint32_t g;
    uint32_t s;
} small_prime;

typedef __int128_t int128_t;

typedef struct {

    int16_t coeffs[COMPRESSED_POOL_SIZE];
    
    int64_t base_case[4];

    int64_t input_stack[STACK_BUF_SIZE];
    int64_t workspace[WORK_BUF_SIZE];
} MemContext;  //19KB

typedef struct {
    MemContext x_ctx;
    MemContext y_ctx;
    int64_t c_coeffs[K_VAL * N_MAX]; 
} MIGD_Output; //50KB

typedef struct {

    int64_t u_hat_num[ANTRAG_D]; 
    int64_t c11[ANTRAG_D];
    int64_t c21[ANTRAG_D];
    int64_t c22[ANTRAG_D];

    MIGD_Output migd_key;

    int valid; 
} PreMatrix_Output; //66KB

typedef struct {
    int64_t coeffs[ANTRAG_D]; 
} poly_int64;

typedef struct {
    int8_t f[ANTRAG_D];
    int8_t g[ANTRAG_D];
    int8_t F[ANTRAG_D];
    int8_t G[ANTRAG_D];
    PreMatrix_Output mat;
} PrivateKey;  //72KB

typedef struct {
    int16_t h[ANTRAG_D];
} PublicKey;

int keygen_fg_impl(secret_key_fg *sk);
int antrag_solve_ntru(int8_t *F, int8_t *G, const int8_t *f, const int8_t *g, unsigned logn, uint32_t *tmp);
int Run_PreMatrix(const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G, PreMatrix_Output *out);

#endif