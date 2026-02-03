#ifndef PREMATRIX_H
#define PREMATRIX_H

#include <stdint.h>
#include "common.h"     
#include "migd.h"       
#include "fft.h"        

// 确保 N_MAX 至少为 512
#ifndef N_MAX
    #include "eigd.h"
#endif

// 如果 eigd.h 中的 N_MAX 小于 ANTRAG_D (512)，我们需要定义自己的常量
// 为了安全，直接定义 PREMATRIX_N 为 512
#define PREMATRIX_N 512

typedef struct {
    // 使用显式的 512 大小，避免依赖 N_MAX 可能带来的隐患
    int64_t u_hat_num[PREMATRIX_N]; 
    int64_t c11[PREMATRIX_N];
    int64_t c21[PREMATRIX_N];
    int64_t c22[PREMATRIX_N];

    MIGD_Output migd_key;

    int valid; 
} PreMatrix_Output;

int Run_PreMatrix(
    const int8_t *f, const int8_t *g, 
    const int8_t *F, const int8_t *G, 
    PreMatrix_Output *out
);

#endif // PREMATRIX_H