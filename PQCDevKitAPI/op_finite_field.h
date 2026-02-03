#ifndef OP_FINITE_FIELD_H
#define OP_FINITE_FIELD_H

#include <stdint.h>
#include "operator_interface.h"

// 对应 switch 的操作码
#define OP_ADD 0
#define OP_SUB 1
#define OP_MUL 2
#define OP_POW 3
#define OP_INV 4

#define MAX_BIGINT_LEN 64 

// -------------------------------------------------------------------------
// 封装接口
// 输入指针 a, b, n 指向 64-bit 端子数据
// len：数据的端子个数 (Limbs)
//      L1 = 4 (对应 32 bytes)
//      L3 = 6 (对应 48 bytes)
//      L5 = 8 (对应 64 bytes)
// -------------------------------------------------------------------------

static inline void ff_add(uint64_t *c, const uint64_t *a, const uint64_t *b, const uint64_t *n, int len) {
    ABORT_IF_FAIL(OP_finite_field(c, OP_ADD, a, b, n, len));
}

static inline void ff_sub(uint64_t *c, const uint64_t *a, const uint64_t *b, const uint64_t *n, int len) {
    ABORT_IF_FAIL(OP_finite_field(c, OP_SUB, a, b, n, len));
}

static inline void ff_mul(uint64_t *c, const uint64_t *a, const uint64_t *b, const uint64_t *n, int len) {
    ABORT_IF_FAIL(OP_finite_field(c, OP_MUL, a, b, n, len));
}

static inline void ff_pow(uint64_t *c, const uint64_t *a, const uint64_t *b, const uint64_t *n, int len) {
    ABORT_IF_FAIL(OP_finite_field(c, OP_POW, a, b, n, len));
}

static inline void ff_inverse(uint64_t *c, const uint64_t *a, const uint64_t *n, int len) {
    ABORT_IF_FAIL(OP_finite_field(c, OP_INV, a, NULL, n, len));
}

#endif // OP_FINITE_FIELD_H