/* tests/test_sign_impl.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "common.h"

// 声明 KeyGen 和 Sign
int crypto_sign_keypair(PublicKey *pk, PrivateKey *sk);
int crypto_sign(uint8_t *sig, size_t *sig_len, const uint8_t *m, size_t mlen, const PrivateKey *sk);

int main() {
    printf("=== Zitaka Signature Test (In-Memory KeyGen) ===\n");

    PublicKey pk;
    PrivateKey sk;

    // 1. 当场生成密钥 (避开文件读写导致的指针/数据丢失问题)
    printf("Generating fresh KeyPair...\n");
    clock_t start = clock();
    if (!crypto_sign_keypair(&pk, &sk)) {
        printf("[FATAL] KeyGen failed.\n");
        return 1;
    }
    printf("KeyGen done in %.4fs.\n", (double)(clock()-start)/CLOCKS_PER_SEC);

    // 2. 准备消息
    const char *msg = "Hello Zitaka!";
    size_t mlen = strlen(msg);
    
    // 3. 执行签名
    size_t sig_len = 2048; 
    uint8_t *sig = malloc(sig_len);
    
    printf("Signing...\n");
    int ret = crypto_sign(sig, &sig_len, (const uint8_t*)msg, mlen, &sk);

    if (ret == 1) {
        printf("[SUCCESS] Signature generated!\n");
        printf("Signature Length: %lu bytes\n", sig_len);
        
        // 验证压缩头
        printf("Header: %02x (Expect 39)\n", sig[40]);
        
        // 简单验证范数 (可选)
        // verify_sig(sig, msg, pk); 
    } else {
        printf("[FAILED] Signature generation failed (ret=%d)\n", ret);
    }

    free(sig);
    return 0;
}