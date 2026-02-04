#include <stdint.h>
#include <stddef.h>
#include "common.h"


size_t compress_sig(void *out, size_t max_out_len, const int16_t *x) {
    uint8_t *buf = (uint8_t *)out;
    size_t n = ANTRAG_D;
    size_t u, v;
    uint32_t acc;
    unsigned acc_len;

    for (u = 0; u < n; u++) {
        if (x[u] < -2047 || x[u] > +2047) {
            return 0;
        }
    }

    acc = 0;
    acc_len = 0;
    v = 0;

    if (buf != NULL) {
        if (v >= max_out_len) return 0;
        buf[v] = 0x30 | 9;
    }
    v++;

    for (u = 0; u < n; u++) {
        int t = x[u];
        unsigned w;

        acc <<= 1;
        if (t < 0) {
            t = -t;
            acc |= 1;
        }
        w = (unsigned)t;

        acc <<= 7;
        acc |= w & 127u;
        w >>= 7;

        acc_len += 8;

        acc <<= (w + 1);
        acc |= 1;
        acc_len += w + 1;

        while (acc_len >= 8) {
            acc_len -= 8;
            if (buf != NULL) {
                if (v >= max_out_len) return 0;
                buf[v] = (uint8_t)(acc >> acc_len);
            }
            v++;
        }
    }

    if (acc_len > 0) {
        if (buf != NULL) {
            if (v >= max_out_len) return 0;
            buf[v] = (uint8_t)(acc << (8 - acc_len));
        }
        v++;
    }

    return v;
}

size_t decompress_sig(int16_t *x, const uint8_t *in, size_t max_in_len) {
    size_t n = ANTRAG_D;
    uint32_t acc;
    unsigned acc_len;
    size_t u, v;

    if (max_in_len < 1) return 0;
    
    // 1. 检查 Header (Falcon-512 = 0x39)
    if (in[0] != (0x30 | 9)) return 0;
    
    acc = 0;
    acc_len = 0;
    v = 1; // 从第1个字节开始读

    for (u = 0; u < n; u++) {
        unsigned int sign;
        unsigned int low;
        uint32_t high;

        // 读 1 bit 符号位
        while (acc_len < 1) {
            if (v >= max_in_len) return 0;
            acc = (acc << 8) | in[v++];
            acc_len += 8;
        }
        sign = (acc >> (acc_len - 1)) & 1;
        acc_len --;

        // 读 7 bits 低位
        while (acc_len < 7) {
            if (v >= max_in_len) return 0;
            acc = (acc << 8) | in[v++];
            acc_len += 8;
        }
        low = (acc >> (acc_len - 7)) & 127;
        acc_len -= 7;

        // 读高位一元编码 (0...01)
        high = 0;
        for (;;) {
            while (acc_len < 1) {
                if (v >= max_in_len) return 0;
                acc = (acc << 8) | in[v++];
                acc_len += 8;
            }
            if ((acc >> (acc_len - 1)) & 1) {
                acc_len --;
                break;
            }
            high ++;
            acc_len --;
            if (high > 31) return 0; // 防止恶意构造导致溢出
        }

        // 组合
        int val = (int)low | ((int)high << 7);
        if (sign) {
            val = -val;
        }
        
        // 范围检查
        if (u == 0 && val == 0 && sign != 0) return 0; // Falcon 规范：-0 是非法的
        x[u] = (int16_t)val;
    }

    return v;
}