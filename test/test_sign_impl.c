#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "common.h"

int crypto_sign_keypair(PublicKey *pk, PrivateKey *sk);
int crypto_sign(uint8_t *sig, size_t *sig_len, const uint8_t *m, size_t mlen, const PrivateKey *sk);
int crypto_verify(const uint8_t *sig, size_t sig_len, const uint8_t *m, size_t mlen, const PublicKey *pk);


int main(void) {
    clock_t start;
    printf("========================================\n");
    printf("   Zitaka / ANTRAG_Z Signature Test     \n");
    printf("========================================\n");

    int ret;
    PublicKey pk;
    PrivateKey sk;
    uint8_t msg[] = "Hello, Zitaka ANTRAG_Z Hardware Implementation!";
    size_t mlen = strlen((char*)msg);
    uint8_t sig[2048];
    size_t sig_len = sizeof(sig);

    // 1. Key Generation
    printf("\n[1] Generating Keypair...\n");
    ret = crypto_sign_keypair(&pk, &sk);

    if (ret != 1) {
        fprintf(stderr, "KeyGen Failed!\n");
        return 1;
    }
    printf("KeyGen Success. PK Hash (preview): %d\n", pk.h[0]);

    // 2. Signing
    printf("\n[2] Signing Message...\n");

    sig_len = sizeof(sig); 
    ret = crypto_sign(sig, &sig_len, msg, mlen, &sk);

    if (ret != 1) {
        fprintf(stderr, "Sign Failed (ret=%d)!\n", ret);
        return 1;
    }
    printf("Sign Success. Signature Length: %zu bytes\n", sig_len);

    // 3. Verification
    printf("\n[3] Verifying Signature...\n");

    ret = crypto_verify(sig, sig_len, msg, mlen, &pk);

    if (ret == 1) {
        printf("\n>>> VERIFICATION SUCCESSFUL <<<\n");
    } else {
        printf("\n>>> VERIFICATION FAILED <<<\n");
        return 1;
    }

    // 4. Corrupt Test (Optional)
    printf("\n[4] Corruption Test...\n");
    sig[50] ^= 0x01;
    ret = crypto_verify(sig, sig_len, msg, mlen, &pk);
    if (ret == 0) {
        printf("Corruption correctly detected.\n");
    } else {
        printf("ERROR: Corrupted signature passed verification!\n");
    }

    return 0;
}