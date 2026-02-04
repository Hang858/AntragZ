#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "fft.h" 
#include "poly.h"
#include "file_utils.h"

int crypto_sign_keypair(PublicKey *pk, PrivateKey *sk) {
    memset(sk, 0, sizeof(PrivateKey));
    memset(pk, 0, sizeof(PublicKey));
    if (!pk || !sk) return 0;
    
    printf("[Zitaka] KeyGen Step 1: Generating f, g...\n");

    secret_key_fg fg_ctx;
    int retries = 500;
    int success = 0;
    int16_t *f_inv_check = malloc(ANTRAG_D * sizeof(int16_t));

    if (!f_inv_check) return 0;

    while(retries--) {
        keygen_fg_impl(&fg_ctx);
        if (poly_inv_mod_q(f_inv_check, fg_ctx.f)) {
            success = 1;
            break;
        }
    }
    free(f_inv_check);

    if (!success) {
        printf("[Error] Failed to generate invertible f within retry limit.\n");
        return 0;
    }
    
    memcpy(sk->f, fg_ctx.f, ANTRAG_D);
    memcpy(sk->g, fg_ctx.g, ANTRAG_D);

    printf("[Zitaka] KeyGen Step 2: Solving NTRU Equation...\n");

    uint32_t *tmp_ntru = malloc(9000 * sizeof(uint32_t)); 
    if (!tmp_ntru) return 0;
    
    if (!antrag_solve_ntru(sk->F, sk->G, sk->f, sk->g, ANTRAG_LOGD, tmp_ntru)) {
        printf("[Error] Failed to solve NTRU equation.\n");
        free(tmp_ntru);
        return 0;
    }
    free(tmp_ntru);

    printf("[Zitaka] KeyGen Step 3: Computing PreMatrix...\n");

    // 3. 计算 PreMatrix (签名所需的分解矩阵)
    if (!Run_PreMatrix(sk->f, sk->g, sk->F, sk->G, &sk->mat)) {
        printf("[Error] PreMatrix computation failed (Cholesky/MIGD error).\n");
        return 0;
    }
    
    printf("[Zitaka] KeyGen Step 4: Computing Public Key h...\n");

    // 4. 生成公钥 h = g * f^-1 mod q
    int16_t *f_inv = malloc(ANTRAG_D * sizeof(int16_t));
    int16_t *g_poly = malloc(ANTRAG_D * sizeof(int16_t));
    
    if (!f_inv || !g_poly) {
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

    printf("[Zitaka] Saving keys to file...\n");
    if (save_to_file("sk.bin", sk, sizeof(PrivateKey))) {
        printf("[Info] Secret Key saved to 'sk.bin' (%zu bytes)\n", sizeof(PrivateKey));
    }
    if (save_to_file("pk.bin", pk, sizeof(PublicKey))) {
        printf("[Info] Public Key saved to 'pk.bin' (%zu bytes)\n", sizeof(PublicKey));
    }

    printf("[Zitaka] KeyPair Generated Successfully.\n");
    return 1;

    printf("[Zitaka] KeyPair Generated Successfully.\n");
    return 1;
}