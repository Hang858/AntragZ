#include <stdlib.h>
#include <string.h>
#include "poly.h"

// 扩展欧几里得求整数逆元
int16_t mod_inverse(int16_t a, int16_t m) {
    int t = 0, newt = 1;
    int r = m, newr = a;
    while (newr != 0) {
        int quotient = r / newr;
        int temp;
        temp = t; t = newt; newt = temp - quotient * newt;
        temp = r; r = newr; newr = temp - quotient * newr;
    }
    if (r > 1) return 0; // 不可逆
    if (t < 0) t += m;
    return t;
}

// 多项式乘法 mod (x^n + 1, q)
// 朴素 O(N^2) 实现
void poly_mul_mod(int16_t *res, const int16_t *a, const int16_t *b) {
    memset(res, 0, ANTRAG_D * sizeof(int16_t));
    for (int i = 0; i < ANTRAG_D; i++) {
        for (int j = 0; j < ANTRAG_D; j++) {
            int k = i + j;
            int32_t val = (int32_t)a[i] * b[j];
            
            // X^N = -1 mod (X^N + 1)
            if (k < ANTRAG_D) {
                res[k] = mod_q(res[k] + val);
            } else {
                res[k - ANTRAG_D] = mod_q(res[k - ANTRAG_D] - val);
            }
        }
    }
}

// 多项式求逆 mod (x^n + 1, q)
// 使用扩展欧几里得算法
int poly_inv_mod_q(int16_t *out, const int8_t *f_in) {
    int n = ANTRAG_D;
    int32_t q = ANTRAG_Q;

    int16_t *u = calloc(n, sizeof(int16_t));
    int16_t *v = calloc(n + 1, sizeof(int16_t));
    int16_t *b = calloc(n, sizeof(int16_t));
    int16_t *c = calloc(n, sizeof(int16_t));
    int16_t *tmp = calloc(n, sizeof(int16_t));

    if (!u || !v || !b || !c || !tmp) goto fail;

    // u = f
    for (int i = 0; i < n; i++) u[i] = mod_q(f_in[i]);
    int deg_u = n - 1;
    while (deg_u >= 0 && u[deg_u] == 0) deg_u--;

    // v = x^n + 1
    memset(v, 0, (n + 1) * sizeof(int16_t));
    v[0] = 1; v[n] = 1; 
    int deg_v = n;

    // b = 1, c = 0
    b[0] = 1; 

    while (deg_u != -1) {
        if (deg_u == 0) {
            int16_t inv = mod_inverse(u[0], q);
            if (inv == 0) goto fail;

            for (int i = 0; i < n; i++) {
                out[i] = mod_q((int32_t)b[i] * inv);
            }
            
            free(u); free(v); free(b); free(c); free(tmp);
            return 1; // 成功
        }

        if (deg_u < deg_v) {
            // swap u, v
            int16_t *ptr; int deg;
            ptr = u; u = v; v = ptr;
            deg = deg_u; deg_u = deg_v; deg_v = deg;
            // swap b, c
            ptr = b; b = c; c = ptr;
        }

        int diff = deg_u - deg_v;
        int scale = mod_q((int32_t)u[deg_u] * mod_inverse(v[deg_v], q));

        // u -= scale * v * x^diff
        for (int i = 0; i <= deg_v; i++) {
            int target_idx = i + diff;
            int32_t sub = mod_q((int32_t)v[i] * scale);
            u[target_idx] = mod_q(u[target_idx] - sub);
        }

        // b -= scale * c * x^diff
        memset(tmp, 0, n * sizeof(int16_t));
        for (int i = 0; i < n; i++) {
            if (c[i] == 0) continue;
            int pos = i + diff;
            int32_t val = mod_q((int32_t)c[i] * scale);
            
            if (pos < n) {
                tmp[pos] = mod_q(tmp[pos] + val);
            } else {
                tmp[pos - n] = mod_q(tmp[pos - n] - val);
            }
        }
        for(int i=0; i<n; i++) {
            b[i] = mod_q(b[i] - tmp[i]);
        }

        while (deg_u >= 0 && u[deg_u] == 0) deg_u--;
    }

fail:
    if(u) free(u); if(v) free(v); if(b) free(b); if(c) free(c); if(tmp) free(tmp);
    return 0;
}