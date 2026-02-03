#include <stdint.h>
#include <string.h>
#include <stddef.h>

// 128-bit 类型支持 (用于乘法高位计算)
typedef unsigned __int128 uint128_t;
// 最大支持64字节大整数运算
#define MAX_BIGINT_LEN 64 

// =========================================================================
//  常数定义 (64-bit 端子)，对应 SQIsign 的三个特定素数
// =========================================================================

// L1: p = 5 * 2^248 - 1 (4 limbs)
static const uint64_t P1_LIMBS[4] = {
    0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 
    0xFFFFFFFFFFFFFFFFULL, 0x04FFFFFFFFFFFFFFULL
};

// L3: p = 65 * 2^376 - 1 (6 limbs)
static const uint64_t P3_LIMBS[6] = {
    0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
    0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0x40FFFFFFFFFFFFFFULL
};

// L5: p = 27 * 2^500 - 1 (8 limbs)
static const uint64_t P5_LIMBS[8] = {
    0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
    0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0x1AFFFFFFFFFFFFFULL
};

// =========================================================================
//  64-bit 核心运算辅助函数
// =========================================================================

static inline void mul64(uint64_t a, uint64_t b, uint64_t *hi, uint64_t *lo) {
    uint128_t res = (uint128_t)a * b;
    *lo = (uint64_t)res;
    *hi = (uint64_t)(res >> 64);
}

static inline uint64_t adc64(uint64_t a, uint64_t b, uint64_t carry, uint64_t *r) {
    uint128_t res = (uint128_t)a + b + carry;
    *r = (uint64_t)res;
    return (uint64_t)(res >> 64);
}

static inline uint64_t add64(uint64_t a, uint64_t b, uint64_t *r) {
    uint64_t res = a + b;
    *r = res;
    return (res < a) ? 1 : 0;
}

static inline uint64_t sbb64(uint64_t a, uint64_t b, uint64_t borrow, uint64_t *r) {
    uint128_t res = (uint128_t)a - b - borrow;
    *r = (uint64_t)res;
    return (uint64_t)(res >> 127);
}

// =========================================================================
//  核心算术逻辑 (按 64-bit Limbs 处理)
// =========================================================================

static uint64_t vec_add(uint64_t *c, const uint64_t *a, const uint64_t *b, int n) {
    uint64_t carry = 0;
    for (int i = 0; i < n; i++) {
        carry = adc64(a[i], b[i], carry, &c[i]);
    }
    return carry;
}

static uint64_t vec_sub(uint64_t *c, const uint64_t *a, const uint64_t *b, int n) {
    uint64_t borrow = 0;
    for (int i = 0; i < n; i++) {
        borrow = sbb64(a[i], b[i], borrow, &c[i]);
    }
    return borrow;
}

static void vec_copy(uint64_t *dst, const uint64_t *src, int n) {
    memcpy(dst, src, n * sizeof(uint64_t));
}

// 模加
static void mod_add(uint64_t *c, const uint64_t *a, const uint64_t *b, const uint64_t *p, int n) {
    uint64_t carry = vec_add(c, a, b, n);
    uint64_t tmp[8];
    uint64_t borrow = vec_sub(tmp, c, p, n);
    if (carry || borrow == 0) {
        vec_copy(c, tmp, n);
    }
}

// 模减
static void mod_sub(uint64_t *c, const uint64_t *a, const uint64_t *b, const uint64_t *p, int n) {
    uint64_t borrow = vec_sub(c, a, b, n);
    if (borrow) {
        vec_add(c, c, p, n);
    }
}

// 蒙哥马利乘法 (CIOS, mu=1)
static void mont_mul(uint64_t *c, const uint64_t *a, const uint64_t *b, const uint64_t *p, int n) {
    uint64_t t[10] = {0}; 
    for (int i = 0; i < n; i++) {
        uint64_t u = 0; 
        for (int j = 0; j < n; j++) {
            uint64_t hi, lo, c1, c2;
            mul64(a[i], b[j], &hi, &lo);
            c1 = add64(t[j], lo, &t[j]);
            c2 = add64(t[j], u, &t[j]);
            u = hi + c1 + c2;
        }
        add64(t[n], u, &t[n]);

        uint64_t m = t[0]; // mu=1
        u = 0;
        uint64_t val_lo, val_hi;
        
        mul64(m, p[0], &val_hi, &val_lo);
        add64(t[0], val_lo, &t[0]);
        u = val_hi + (t[0] < val_lo);

        for (int j = 1; j < n; j++) {
            uint64_t c1, c2;
            mul64(m, p[j], &val_hi, &val_lo);
            c1 = add64(t[j], val_lo, &t[j]);
            c2 = add64(t[j], u, &t[j]);
            u = val_hi + c1 + c2;
        }
        
        uint64_t c1;
        c1 = add64(t[n], u, &t[n]);
        add64(t[n+1], c1, &t[n+1]);

        memmove(t, &t[1], (n + 1) * sizeof(uint64_t));
        t[n+1] = 0;
    }

    uint64_t borrow = vec_sub(c, t, p, n);
    if (borrow) {
        vec_copy(c, t, n);
    }
}

// 模幂
static void mont_pow(uint64_t *c, const uint64_t *a, const uint64_t *exp, const uint64_t *p, int n) {
    int bit_len = n * 64;
    int msb = bit_len - 1;
    while (msb >= 0 && !((exp[msb / 64] >> (msb % 64)) & 1)) {
        msb--;
    }
    if (msb < 0) {
        memset(c, 0, n * sizeof(uint64_t)); 
        return; 
    }
    uint64_t res[8];
    vec_copy(res, a, n);
    for (int i = msb - 1; i >= 0; i--) {
        mont_mul(res, res, res, p, n);
        int limb_idx = i / 64;
        int bit_idx  = i % 64;
        if ((exp[limb_idx] >> bit_idx) & 1) {
            mont_mul(res, res, a, p, n);
        }
    }
    vec_copy(c, res, n);
}

// 模逆: a^(p-2)
static void mont_inv(uint64_t *c, const uint64_t *a, const uint64_t *p, int n) {
    uint64_t exp[8];
    vec_copy(exp, p, n);
    exp[0] -= 2;
    mont_pow(c, a, exp, p, n);
}

// =========================================================================
//  主接口实现
// =========================================================================

int OP_finite_field(uint64_t *c_out, uint8_t opr, const uint64_t *a_in, const uint64_t *b_in, const uint64_t *n_in, int len)
{
    // 1. 安全检查
    if (!c_out || !a_in || !n_in) return -1;

    // 2. 长度校验
    if (len != 4 && len != 6 && len != 8) return -1;
    
    // 3. 选择素数
    const uint64_t *p_limbs;
    if (len == 4) p_limbs = P1_LIMBS;
    else if (len == 6) p_limbs = P3_LIMBS;
    else p_limbs = P5_LIMBS;

    // 4. 本地缓冲区 (避免指针重叠问题)
    uint64_t c[8]; 

    // 5. 执行运算 (直接使用 len 作为循环次数)
    switch (opr) {
        case 0: // ADD
            mod_add(c, a_in, b_in, p_limbs, len);
            break;
        case 1: // SUB
            mod_sub(c, a_in, b_in, p_limbs, len);
            break;
        case 2: // MUL
            mont_mul(c, a_in, b_in, p_limbs, len);
            break;
        case 3: // POW
            mont_pow(c, a_in, b_in, p_limbs, len);
            break;
        case 4: // INV
            mont_inv(c, a_in, p_limbs, len);
            break;
        default:
            return -1;
    }

    // 6. 输出结果
    memcpy(c_out, c, len * sizeof(uint64_t));

    return 0;
}