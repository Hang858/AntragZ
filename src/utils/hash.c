#include <stdint.h>
#include <stddef.h>
#include "common.h"
#include "operator_interface.h"
#include "hash.h"

#define FALCON_Q 12289

void hash_to_point(int16_t *x, const uint8_t *r, size_t r_len, const uint8_t *msg, size_t msg_len) {

    shake256incctx sc; 
    // shake256_inc_init(&sc);
    // shake256_inc_absorb(&sc, r, r_len);
    // shake256_inc_absorb(&sc, msg, msg_len);
    // shake256_inc_finalize(&sc);
    OP_hash_init(OP_ALG_SHAKE256, (void*)sc.ctx, 208);
    OP_hash_absorb(OP_ALG_SHAKE256, (void*)sc.ctx, 208, (void*)r, r_len);
    OP_hash_absorb(OP_ALG_SHAKE256, (void*)sc.ctx, 208, (void*)msg, msg_len);

    size_t n = ANTRAG_D; // 512
    
    while (n > 0) {
        uint8_t buf[2];
        uint32_t w;

        OP_hash_squeeze(OP_ALG_SHAKE256, (void*)sc.ctx, 208, buf, 2);

        w = ((unsigned)buf[0] << 8) | (unsigned)buf[1];
        
        if (w < 61445) {
            while (w >= ANTRAG_Q) {
                w -= ANTRAG_Q;
            }
            *x++ = (int16_t)w;
            n--;
        }
    }
}