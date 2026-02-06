/* test/test_dump_samples.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "common.h"
#include "hash.h"
#include "poly.h"
#include "encode.h"
#include "rng.h"

#define NUM_SAMPLES 10000
int crypto_sign_keypair(PublicKey *pk, PrivateKey *sk);
int crypto_sign(uint8_t *sig, size_t *sig_len, const uint8_t *m, size_t mlen, const PrivateKey *sk);
// s2 = s1 * h + c (mod q)
static void reconstruct_s2(int16_t *s2, const int16_t *s1, const int16_t *c_hash, const PublicKey *pk) {
    int16_t tmp_s2[ANTRAG_D];
    
    // 计算 s1 * h
    poly_mul_mod(tmp_s2, s1, pk->h);

    // s2 = s1 * h + c (mod q, centered)
    for (int i = 0; i < ANTRAG_D; i++) {
        int32_t val = (int32_t)tmp_s2[i] + c_hash[i];
        
        // Center mod q to [-q/2, q/2]
        int32_t r = val % ANTRAG_Q;
        if (r > ANTRAG_Q / 2) r -= ANTRAG_Q;
        if (r < -ANTRAG_Q / 2) r += ANTRAG_Q;
        
        s2[i] = (int16_t)r;
    }
}

int main(void) {
    printf("=== Generating %d Samples for Covariance Check ===\n", NUM_SAMPLES);

    // 1. 生成密钥对
    PublicKey pk;
    PrivateKey sk;
    if (!crypto_sign_keypair(&pk, &sk)) {
        fprintf(stderr, "Keygen failed!\n");
        return 1;
    }

    // 2. 打开输出文件
    FILE *fp = fopen("samples.bin", "wb");
    if (!fp) {
        perror("Failed to open samples.bin");
        return 1;
    }

    uint8_t msg[] = "TestMessage";
    uint8_t sig[2048];
    size_t sig_len;
    
    int16_t s1[ANTRAG_D];
    int16_t s2[ANTRAG_D];
    int16_t c_hash[ANTRAG_D];
    
    // 缓冲区：存储 [s1, s2] 拼接向量
    // 每个样本大小 = 512 * 2 * 2 bytes = 2048 bytes
    int16_t sample_vec[2 * ANTRAG_D];

    for (int i = 0; i < NUM_SAMPLES; i++) {
        sig_len = sizeof(sig);
        
        // 3. 签名
        // 为了统计随机性，每次签名都使用新的随机数（crypto_sign 内部会生成）
        if (crypto_sign(sig, &sig_len, msg, sizeof(msg), &sk) != 1) {
            fprintf(stderr, "Sign failed at iter %d\n", i);
            continue;
        }

        // 4. 解析签名数据
        // 4.1 恢复 Hash (c)
        hash_to_point(c_hash, sig, 40, msg, sizeof(msg));
        
        // 4.2 解压 s1
        if (decompress_sig(s1, sig + 40, sig_len - 40) == 0) {
            fprintf(stderr, "Decompress failed at iter %d\n", i);
            continue;
        }

        // 4.3 重构 s2
        reconstruct_s2(s2, s1, c_hash, &pk);

        // 5. 拼接并写入文件 [s1_0...s1_n, s2_0...s2_n]
        memcpy(sample_vec, s1, ANTRAG_D * sizeof(int16_t));
        memcpy(sample_vec + ANTRAG_D, s2, ANTRAG_D * sizeof(int16_t));
        
        fwrite(sample_vec, sizeof(int16_t), 2 * ANTRAG_D, fp);

        if (i % 1000 == 0) {
            printf("Generated %d/%d samples...\r", i, NUM_SAMPLES);
            fflush(stdout);
        }
    }

    printf("\nDone. Data saved to 'samples.bin'.\n");
    fclose(fp);
    return 0;
}