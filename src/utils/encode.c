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
    
    if (in[0] != (0x30 | 9)) return 0;
    
    acc = 0;
    acc_len = 0;
    v = 1;

    for (u = 0; u < n; u++) {
        unsigned int sign;
        unsigned int low;
        uint32_t high;

        while (acc_len < 1) {
            if (v >= max_in_len) return 0;
            acc = (acc << 8) | in[v++];
            acc_len += 8;
        }
        sign = (acc >> (acc_len - 1)) & 1;
        acc_len --;

        while (acc_len < 7) {
            if (v >= max_in_len) return 0;
            acc = (acc << 8) | in[v++];
            acc_len += 8;
        }
        low = (acc >> (acc_len - 7)) & 127;
        acc_len -= 7;

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
            if (high > 31) return 0;
        }

        int val = (int)low | ((int)high << 7);
        if (sign) {
            val = -val;
        }
        
        if (u == 0 && val == 0 && sign != 0) return 0;
        x[u] = (int16_t)val;
    }

    return v;
}

void encode_public_key(uint8_t *out, const int16_t *h) {

    out[0] = (uint8_t)ANTRAG_LOGD; 
    uint8_t *buf = out + 1;
    
    for (int i = 0; i < ANTRAG_D; i += 4) {
        uint16_t t0 = h[i+0] & 0x3FFF;
        uint16_t t1 = h[i+1] & 0x3FFF;
        uint16_t t2 = h[i+2] & 0x3FFF;
        uint16_t t3 = h[i+3] & 0x3FFF;

        buf[0] = (uint8_t)(t0);
        buf[1] = (uint8_t)(t0 >> 8) | (uint8_t)(t1 << 6);
        buf[2] = (uint8_t)(t1 >> 2);
        buf[3] = (uint8_t)(t1 >> 10) | (uint8_t)(t2 << 4);
        buf[4] = (uint8_t)(t2 >> 4);
        buf[5] = (uint8_t)(t2 >> 12) | (uint8_t)(t3 << 2);
        buf[6] = (uint8_t)(t3 >> 6);
        
        buf += 7;
    }
}


int decode_public_key(int16_t *h, const uint8_t *in, size_t len) {
    if (len < 897) return 0;
    
    if (in[0] != ANTRAG_LOGD) return 0; 
    
    const uint8_t *buf = in + 1;
    for (int i = 0; i < ANTRAG_D; i += 4) {
        uint8_t b0 = buf[0];
        uint8_t b1 = buf[1];
        uint8_t b2 = buf[2];
        uint8_t b3 = buf[3];
        uint8_t b4 = buf[4];
        uint8_t b5 = buf[5];
        uint8_t b6 = buf[6];
        buf += 7;

        h[i+0] = (int16_t)(b0 | ((b1 & 0x3F) << 8));
        h[i+1] = (int16_t)((b1 >> 6) | (b2 << 2) | ((b3 & 0x0F) << 10));
        h[i+2] = (int16_t)((b3 >> 4) | (b4 << 4) | ((b5 & 0x03) << 12));
        h[i+3] = (int16_t)((b5 >> 2) | (b6 << 6));
    }
    return 1;
}