#ifndef MATH_NTT_H
#define MATH_NTT_H

#include <stdint.h>
#include <stddef.h>
#include "modp.h"

// 预计算生成元表
void modp_mkgm2(uint32_t * gm, uint32_t * igm, unsigned logn,
    uint32_t g, uint32_t p, uint32_t p0i);

// NTT 变换
void modp_NTT2_ext(uint32_t *a, size_t stride, const uint32_t *gm, unsigned logn,
    uint32_t p, uint32_t p0i);

// 逆 NTT 变换
void modp_iNTT2_ext(uint32_t *a, size_t stride, const uint32_t *igm, unsigned logn,
    uint32_t p, uint32_t p0i);

// 简化宏
#define modp_NTT2(a, gm, logn, p, p0i)   modp_NTT2_ext(a, 1, gm, logn, p, p0i)
#define modp_iNTT2(a, igm, logn, p, p0i) modp_iNTT2_ext(a, 1, igm, logn, p, p0i)

#endif