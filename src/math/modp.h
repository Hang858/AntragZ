#ifndef MODP_H
#define MODP_H

#include <stdint.h>
#include <stddef.h>

// 设置 int32_t 值到 mod p 域
static inline uint32_t modp_set(int32_t x, uint32_t p) {
    uint32_t w = (uint32_t)x;
    w += p & -(w >> 31);
    return w;
}

// 将 mod p 域的值规范化到 [-p/2, p/2]
static inline int32_t modp_norm(uint32_t x, uint32_t p) {
    return (int32_t)(x - (p & (((x - ((p + 1) >> 1)) >> 31) - 1)));
}

static inline uint32_t modp_add(uint32_t a, uint32_t b, uint32_t p) {
    uint32_t d = a + b - p;
    d += p & -(d >> 31);
    return d;
}

static inline uint32_t modp_sub(uint32_t a, uint32_t b, uint32_t p) {
    uint32_t d = a - b;
    d += p & -(d >> 31);
    return d;
}

// 蒙哥马利模乘
static inline uint32_t modp_montymul(uint32_t a, uint32_t b, uint32_t p, uint32_t p0i) {
    uint64_t z = (uint64_t)a * (uint64_t)b;
    uint64_t w = ((z * p0i) & (uint64_t)0x7FFFFFFF) * p;
    uint32_t d = (uint32_t)((z + w) >> 31) - p;
    d += p & -(d >> 31);
    return d;
}

static inline uint32_t modp_R(uint32_t p) { return ((uint32_t)1 << 31) - p; }

// 计算 -1/p mod 2^31
static inline uint32_t modp_ninv31(uint32_t p) {
    uint32_t y = 2 - p;
    for(int i=0; i<4; i++) y *= 2 - p * y;
    return (uint32_t)0x7FFFFFFF & -y;
}

// 计算 R^2 mod p
static inline uint32_t modp_R2(uint32_t p, uint32_t p0i) {
    uint32_t z = modp_R(p);
    z = modp_add(z, z, p); // 2*R
    for(int i=0; i<5; i++) z = modp_montymul(z, z, p, p0i); // ^32
    return (z + (p & -(z & 1))) >> 1; // 2^62
}

// 计算 R^x mod p
static inline uint32_t modp_Rx(unsigned x, uint32_t p, uint32_t p0i, uint32_t R2) {
    int i;
    uint32_t r, z;
    x --;
    r = R2;
    z = modp_R(p);
    for (i = 0; (1U << i) <= x; i ++) {
        if ((x & (1U << i)) != 0) {
            z = modp_montymul(z, r, p, p0i);
        }
        r = modp_montymul(r, r, p, p0i);
    }
    return z;
}

// 模除 a/b mod p
static inline uint32_t modp_div(uint32_t a, uint32_t b, uint32_t p, uint32_t p0i, uint32_t R) {
    uint32_t z, e;
    int i;
    e = p - 2;
    z = R;
    for (i = 30; i >= 0; i --) {
        uint32_t z2;
        z = modp_montymul(z, z, p, p0i);
        z2 = modp_montymul(z, b, p, p0i);
        z ^= (z ^ z2) & -(uint32_t)((e >> i) & 1);
    }
    z = modp_montymul(z, 1, p, p0i);
    return modp_montymul(a, z, p, p0i);
}

#endif // MODP_H