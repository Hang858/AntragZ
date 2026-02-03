#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include "operator_interface.h" 
#include "sm3.h"
#include "sha256.h"
#include "sha3.h"

int OP_hash_init(uint8_t alg, void *s, int s_len) {
    switch (alg)
    {
    case OP_ALG_SHA256:
        if ((size_t)s_len < sizeof(SHA256_CTX)) {
            return OP_FAILURE;
        }
        SHA256_Init((SHA256_CTX *)s);
        return OP_SUCCESS;
    case OP_ALG_SM3:
        if ((size_t)s_len < sizeof(SM3_STATE)) {
            return OP_FAILURE;
        }
        SM3_init((SM3_STATE*)s);
        return OP_SUCCESS;
    case OP_ALG_SHAKE128:
    case OP_ALG_SHAKE256:
    case OP_ALG_SHA3_256:
    case OP_ALG_SHA3_384:
    case OP_ALG_SHA3_512:
        if (s_len < OP_STATE_SIZE_SHA3) {
            return OP_FAILURE;
        }
        memset(s, 0, s_len);
        break;
    default:
        return OP_FAILURE;
    }
    return 0;
}

int OP_hash_absorb(uint8_t alg, void *s, int s_len, const void *input, int input_len) {
    (void)s_len;
    if (alg == OP_ALG_SHA256) {
        SHA256_CTX *state = (SHA256_CTX *)s;
        SHA256_Update(state, input, input_len);
        return 0;
    } else if (alg == OP_ALG_SM3) {
        SM3_STATE *state = (SM3_STATE *)s;
        SM3_process(state, (unsigned char *)input, input_len);
        return 0;
    } else {
        uint64_t *state = (uint64_t *)s;
        unsigned int rate;
        uint8_t padding;

        if (get_alg_params(alg, &rate, &padding) != 0) return -1;
        
        unsigned int pos = GET_POS(state);
        pos = oph_absorb(state, pos, rate, (const uint8_t*)input, input_len);
        SET_POS(state, pos);
        return 0;
    }
}

int OP_hash_squeeze(uint8_t alg, void *s, int s_len, void *output, int output_len) {
    (void)s_len;
    if (alg == OP_ALG_SHA256) {
        SHA256_CTX *state = (SHA256_CTX *)s;
        SHA256_Final((unsigned char *)output, state);
        return 0;
    } else if (alg == OP_ALG_SM3) {
        SM3_STATE *state = (SM3_STATE *)s;
        SM3_done(state, (unsigned char *)output);
        return 0;
    } else {
        uint64_t *state = (uint64_t *)s;
        uint8_t *out = (uint8_t *)output;

        unsigned int rate;
        uint8_t padding;
        
        if (get_alg_params(alg, &rate, &padding) != 0) return -1;

        unsigned int pos = GET_POS(state);
        if (!IS_FINALIZED(state)) {
            ((uint8_t *)state)[pos] ^= padding;
            ((uint8_t *)state)[rate - 1] ^= 0x80;
            
            oph_permute(state);
            
            pos = 0;
            SET_FINALIZED(state);
        }
        pos = oph_squeeze(out, output_len, state, pos, rate);
        SET_POS(state, pos);
        return 0;
    }
}

/**
 * @brief OP_hash 软件参考
 * 
 */
int OP_hash(uint8_t alg, uint8_t mode, int n, const void *input, int input_len, uint8_t link_count, void *output) {
    // 基础参数检查
    if (input == NULL || output == NULL || input_len < 0 || n <= 0) {
        return -1;
    }
    if (mode > 0) {
        if (input_len<34) {
            return -1;
        }
        // 哈希链模式
        uint8_t *ptr=(uint8_t *)input;
        uint8_t *hash=ptr+input_len-32;
        uint8_t *addr=hash-1;
        for (uint8_t i = 0; i < link_count; i++) {
            if (OP_hash(alg, 0, n, input, input_len, 0, hash) != 0) {
                return -1;
            }
            *addr = *addr+1;
        }
        // 最终输出结果
        memcpy(output, hash, 32);
        return 0;
    }
    if (alg == OP_ALG_SHA256) {
        // 调用 SHA-256 实现
        if (n != 32) {
            return -1; // SHA-256 输出长度固定为32字节
        }
        SHA256_CTX ctx;
        SHA256_Init(&ctx);
        SHA256_Update(&ctx, input, input_len);
        SHA256_Final(output, &ctx);
        return 0;       
    } else if (alg == OP_ALG_SM3) {
        // 调用 SM3 实现
        if (n != 32) {
            return -1; // SM3 输出长度固定为32字节
        }
        SM3_256((unsigned char *)input, input_len, (unsigned char *)output);
        return 0;
    }

    uint64_t s[25];
    uint8_t *out_ptr = (uint8_t *)output;
    const uint8_t *in_ptr = (const uint8_t *)input;
    unsigned int rate;
    uint8_t padding;

    // 1. 根据算法选择参数
    switch (alg) {
        case OP_ALG_SHAKE128:
            rate = OPH_SHAKE128_RATE;
            padding = 0x1F;
            break;
        case OP_ALG_SHAKE256:
            rate = OPH_SHAKE256_RATE;
            padding = 0x1F;
            break;
        case OP_ALG_SHA3_256:
            if (n != 32) return -1;
            rate = OPH_SHA3_256_RATE;
            padding = 0x06;
            break;
        case OP_ALG_SHA3_384:
            if (n != 48) return -1;
            rate = OPH_SHA3_384_RATE;
            padding = 0x06;
            break;
        case OP_ALG_SHA3_512:
            if (n != 64) return -1;
            rate = OPH_SHA3_512_RATE;
            padding = 0x06;
            break;
        default:
            return -1; // 不支持的算法
    }

    // 2. 核心逻辑执行
    
    if (alg == OP_ALG_SHA3_256 || alg == OP_ALG_SHA3_512 || alg == OP_ALG_SHA3_384) {
        oph_absorb_once(s, rate, in_ptr, (size_t)input_len, padding);
        oph_permute(s);
        
        unsigned int i;
        for (i = 0; i < (unsigned int)n / 8; i++) {
            oph_store64(out_ptr + 8 * i, s[i]);
        }
        return 0;
    }

    oph_init(s);
    unsigned int pos = 0;
    pos = oph_absorb(s, pos, rate, in_ptr, (size_t)input_len);
    oph_finalize(s, pos, rate, padding);
    
    pos = rate; 
    oph_squeeze(out_ptr, (size_t)n, s, pos, rate);

    return 0;
}
