#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "operator_interface.h"

// ============================================================================
// 参数定义区: Q = 268435457
// ============================================================================

#define BIG_Q         268435457   
#define BIG_INV_Q     4026531839  
// N=256 的逆元保持标准型，因为 OP_intt256 最后是: res * inv_n * R^-1
// 我们希望结果是 Standard，所以输入是 Mont(res), inv_n 必须是 Standard
// 1/256 mod Q = 267386881
#define BIG_INV_N_256 267386881   

// 1/2 的 Montgomery 形式 (用于最后的 ntt_512_inverse_wrapper)
// 1/2 * R mod Q = 134217729 * 16 mod Q = 2147483664
#define BIG_INV_2_MONT 2147483664 

// R^2 mod Q (用于输入转 Montgomery)
// R = 16, R^2 = 256
#define R_SQUARED_MOD_Q 256


// 【重要】请将 Python 脚本生成的表粘贴到这里
// 必须是带 (* R) 的版本
static const int32_t BIG_OMEGA_TABLE_256[] = { /* ... 粘贴脚本输出 ... */ };
static const int32_t PHI_512_FWD[256] = { /* ... 粘贴脚本输出 ... */ };
static const int32_t PHI_512_INV[256] = { /* ... 粘贴脚本输出 ... */ };

// ----------------------------------------------------------------------------
// 辅助函数
// ----------------------------------------------------------------------------
static inline int32_t LocalMontReduce(int64_t a) {
    int64_t t = (int64_t)((int32_t)a * BIG_INV_Q);
    t = (a - t * BIG_Q) >> 32;
    return (int32_t)(t < 0 ? t + BIG_Q : t);
}

// 512 前向变换 (输入 Mont -> 输出 Mont)
static void ntt_512_forward_wrapper(int32_t *buf) {
    int32_t even[256], odd[256];
    for (int i = 0; i < 256; i++) {
        even[i] = buf[2 * i];
        odd[i]  = buf[2 * i + 1];
    }

    // 算子内部全是 MontMul，保持 Mont 域
    OP_ntt256(even, even, 0); 
    OP_ntt256(odd, odd, 0);

    for (int i = 0; i < 256; i++) {
        // PHI 表现在是 Mont 形式
        // odd[i] (Mont) * PHI (Mont) * R^-1 = new_term (Mont)
        int32_t t = LocalMontReduce((int64_t)odd[i] * PHI_512_FWD[i]);
        int32_t u = even[i];

        int32_t val1 = u + t;
        if (val1 >= BIG_Q) val1 -= BIG_Q;
        buf[i] = val1;

        int32_t val2 = u - t;
        if (val2 < 0) val2 += BIG_Q;
        buf[i + 256] = val2;
    }
}

// 512 逆变换 (输入 Mont -> 输出 Standard)
static void ntt_512_inverse_wrapper(int32_t *buf) {
    int32_t even[256], odd[256];

    for (int i = 0; i < 256; i++) {
        int32_t u = buf[i];
        int32_t v = buf[i + 256];

        int32_t sum = u + v;
        if (sum >= BIG_Q) sum -= BIG_Q;
        even[i] = sum; 

        int32_t diff = u - v;
        if (diff < 0) diff += BIG_Q;
        // diff (Mont) * PHI_INV (Mont) * R^-1 = odd (Mont)
        odd[i] = LocalMontReduce((int64_t)diff * PHI_512_INV[i]);
    }

    // 调用硬件算子
    // 注意：OP_intt256 内部最后一步是 * inv_n
    // 我们传入 Standard 的 inv_n，所以输出 = Val(Mont) * inv_n(Std) * R^-1 = Val(Std) * inv_n
    // 即输出变成了 Standard 域！
    OP_intt256(even, even, 1);
    OP_intt256(odd, odd, 1);

    // Merge (此时 even, odd 已经是 Standard 了)
    for (int i = 0; i < 256; i++) {
        buf[2 * i]     = even[i];
        buf[2 * i + 1] = odd[i];
    }
    
    // 最后的 1/2 缩放
    // 输入 buf 是 Standard
    // 我们希望结果是 Standard
    // 技巧：LocalMontReduce(Std * Mont(1/2)) = Std * (1/2 * R) * R^-1 = Std/2
    for (int i = 0; i < 512; i++) {
        buf[i] = LocalMontReduce((int64_t)buf[i] * BIG_INV_2_MONT);
    }
}

void poly_mul_int8_int16_to_32_acc_ntt(int32_t *res, const int8_t *a, const int16_t *b, int n) {
    if (n != 512) return;

    // Init: inv_n 必须传 Standard 形式，因为 OP_ntt256 内部逻辑是 val * inv_n
    OP_ntt256_init(BIG_OMEGA_TABLE_256, 2, BIG_Q, BIG_INV_Q, BIG_INV_N_256, 0);
    
    int32_t *buf_a = malloc(512 * sizeof(int32_t));
    int32_t *buf_b = malloc(512 * sizeof(int32_t));
    
    for (int i = 0; i < n; i++) {
        int32_t val_a = (int32_t)a[i];
        if (val_a < 0) val_a += BIG_Q;
        // 转 Montgomery: val * R = MontReduce(val * R^2)
        // 这一步非常依赖 LocalMontReduce 的正确性
        buf_a[i] = LocalMontReduce((int64_t)val_a * R_SQUARED_MOD_Q);

        int32_t val_b = (int32_t)b[i];
        if (val_b < 0) val_b += BIG_Q;
        buf_b[i] = LocalMontReduce((int64_t)val_b * R_SQUARED_MOD_Q);
    }

    ntt_512_forward_wrapper(buf_a);
    ntt_512_forward_wrapper(buf_b);

    OP_cwm(buf_a, buf_a, buf_b, 256, BIG_Q, 0); 
    OP_cwm(buf_a + 256, buf_a + 256, buf_b + 256, 256, BIG_Q, 0); 

    ntt_512_inverse_wrapper(buf_a); 

    for (int i = 0; i < n; i++) {
        int32_t val = buf_a[i];
        // 中心化还原
        if (val > (BIG_Q >> 1)) val -= BIG_Q;
        res[i] += val;
    }

    free(buf_a);
    free(buf_b);
}