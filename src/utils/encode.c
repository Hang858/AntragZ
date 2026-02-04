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