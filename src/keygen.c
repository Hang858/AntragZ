#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "fft.h" 
#include "poly.h"
#include "file_utils.h"
#include "encode.h"

int crypto_sign_keypair(PublicKey *pk, PrivateKey *sk) {
    // 初始化密钥
    memset(sk, 0, sizeof(PrivateKey));
    memset(pk, 0, sizeof(PublicKey));
    if (!pk || !sk) return 0;


    secret_key_fg fg_ctx;
    int retries = 500;
    int success = 0;
    int16_t *f_inv_check = malloc(ANTRAG_D * sizeof(int16_t));

    if (!f_inv_check) {
        ZITAKA_LOG("Memory allocation failed: f_inv_check");
        return 0;
    }

    uint32_t *tmp_ntru = malloc(9000 * sizeof(uint32_t)); 
    if (!tmp_ntru) {
        ZITAKA_LOG("Memory allocation failed: tmp_ntru");
        free(f_inv_check);
        return 0;
    }
    ZITAKA_LOG("Starting key generation loop...");
    // 重试机制
    while(retries--) {
        // 调用 VectorGen 生成 f, g
        keygen_fg_impl(&fg_ctx);
        // 检查 f 是否可逆
        if (!poly_inv_mod_q(f_inv_check, fg_ctx.f)) {
            continue;
        }
        // 求解 NTRUSolve 生成 F,G
        if (antrag_solve_ntru(sk->F, sk->G, fg_ctx.f, fg_ctx.g, ANTRAG_LOGD, tmp_ntru)) {
            success = 1;
            ZITAKA_LOG("NTRU solve successful at retry %d", 500 - retries);
            break;
        }
    }

    free(f_inv_check);
    free(tmp_ntru);

    if (!success) {
        ZITAKA_LOG("Keygen failed after max retries.");
        return 0;
    }
    // 填充私钥中的 f, g
    memcpy(sk->f, fg_ctx.f, ANTRAG_D);
    memcpy(sk->g, fg_ctx.g, ANTRAG_D);

    // PreMatrix 生成 u_hat 和 A，存入 sk->mat
    if (!Run_PreMatrix(sk->f, sk->g, sk->F, sk->G, &sk->mat)) {
        ZITAKA_LOG("PreMatrix generation failed.");
        return 0;
    }
    
    // 生成公钥 h = g * f^-1 mod q
    int16_t *f_inv = malloc(ANTRAG_D * sizeof(int16_t));
    int16_t *g_poly = malloc(ANTRAG_D * sizeof(int16_t));

    if (!f_inv || !g_poly) {
        ZITAKA_LOG("Memory allocation failed for PK generation.");
        free(f_inv); free(g_poly);
        return 0;
    }

    poly_inv_mod_q(f_inv, sk->f);
    for(int i=0; i<ANTRAG_D; i++) {
        g_poly[i] = mod_q(sk->g[i]);
    }

    poly_mul_mod(pk->h, f_inv, g_poly);

    free(f_inv);
    free(g_poly);
    uint8_t pk_bytes[PK_BYTES];
    encode_public_key(pk_bytes, pk->h);
    ZITAKA_LOG("Saving keys to file...");

    // 公钥私钥保存在本地
#ifdef ZITAKA_DEBUG
    if (save_to_file("sk.bin", sk, sizeof(PrivateKey))) {
        ZITAKA_LOG("Secret Key saved to 'sk.bin' (%zu bytes)", sizeof(PrivateKey));
    }
    if (save_to_file("pk.bin", pk_bytes, PK_BYTES)) {
        ZITAKA_LOG("Public Key saved to 'pk.bin' (%zu bytes)", PK_BYTES);
    }
#endif

    ZITAKA_LOG("KeyPair Generated Successfully.");
    return 1;
}
