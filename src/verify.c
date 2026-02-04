/* src/sign/verify.c */
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "utils/hash.h"
#include <stdlib.h>

// 声明外部工具函数
size_t decompress_sig(int16_t *x, const uint8_t *in, size_t max_in_len);

// 将 mod q 的结果规约到 [-q/2, q/2]
static inline int16_t center_mod(int32_t a) {
    int32_t q = ANTRAG_Q;
    int32_t r = a % q;
    if (r > q / 2) r -= q;
    if (r < -q / 2) r += q;
    return (int16_t)r;
}

int crypto_verify(const uint8_t *sig, size_t sig_len, const uint8_t *m, size_t mlen, const PublicKey *pk) {
    int n = ANTRAG_D;
    int16_t s1[ANTRAG_D];
    int16_t c_hash[ANTRAG_D];
    int16_t s2[ANTRAG_D];

    // 1. 检查签名基本长度 (Salt 40 字节 + 至少 1 字节 Header)
    if (sig_len < 41) return 0;

    // 2. 重新计算挑战值 c = H(r || msg)
    // Salt 存储在签名的前 40 字节
    hash_to_point(c_hash, sig, 40, m, mlen);

    // 3. 解压缩 s1
    size_t bytes_read = decompress_sig(s1, sig + 40, sig_len - 40);
    if (bytes_read == 0) return 0; // 解压缩失败

    // 4. 计算 s2 = s1 * h + c (mod q)
    // 注意：这里需要进行多项式卷积乘法
    int32_t *tmp_s2 = calloc(n, sizeof(int32_t));
    if (!tmp_s2) return 0;

    for (int i = 0; i < n; i++) {
        if (s1[i] == 0) continue;
        for (int j = 0; j < n; j++) {
            int k = i + j;
            int32_t val = (int32_t)s1[i] * pk->h[j];
            if (k < n) {
                tmp_s2[k] += val;
            } else {
                tmp_s2[k - n] -= val; // X^n = -1
            }
        }
    }

    // 5. 计算最终范数并检查
    // 按照算法描述：s2 = (s1*h + c) mod q
    double norm_sq = 0;
    double beta_sq = 34715664.0; // 设计文档中的参数 β^2

    for (int i = 0; i < n; i++) {
        // 计算 s2[i] = (tmp_s2[i] + c_hash[i]) mod q，并规约到中心
        s2[i] = center_mod(tmp_s2[i] + c_hash[i]);
        
        // 累加范数
        norm_sq += (double)s1[i] * s1[i];
        norm_sq += (double)s2[i] * s2[i];
    }

    free(tmp_s2);

    // [DEBUG] 打印验证端的范数
    // printf("[VERIFY] NormSq: %.0f, Limit: %.0f\n", norm_sq, beta_sq);

    if (norm_sq <= beta_sq) {
        return 1; // 验证通过
    } else {
        return 0; // 范数超限
    }
}