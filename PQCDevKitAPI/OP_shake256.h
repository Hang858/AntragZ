#if !defined(SHAKE256_H_)
#define SHAKE256_H_

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include "operator_interface.h"

static inline void shake256_hash(void* input, size_t input_len, uint8_t* output, size_t output_len){
    ABORT_IF_FAIL(OP_hash(OP_ALG_SHAKE256, 0, output_len, input, input_len, 0, output));
}

#define OP_STATE_SIZE_SHA3 208
typedef struct
{
    uint8_t state[OP_STATE_SIZE_SHA3];  
} SHAKE256_STATE;

static inline void shake256_init(SHAKE256_STATE* state){
    ABORT_IF_FAIL(OP_hash_init(OP_ALG_SHAKE256, state, sizeof(SHAKE256_STATE)));
}

static inline void shake256_absorb(SHAKE256_STATE* state, unsigned char* buf, int len){
    ABORT_IF_FAIL(OP_hash_absorb(OP_ALG_SHAKE256, state, sizeof(SHAKE256_STATE), buf, len));
}

static inline void shake256_squeeze(SHAKE256_STATE* state, unsigned char* output, size_t output_len){
    ABORT_IF_FAIL(OP_hash_squeeze(OP_ALG_SHAKE256, state, sizeof(SHAKE256_STATE), output, output_len));
}

static inline void shake256_hash_chain(void* input, size_t input_len, uint8_t* output, size_t output_len, uint8_t link_count){
    ABORT_IF_FAIL(OP_hash(OP_ALG_SHAKE256, 2, output_len, input, input_len, link_count, output));
}

#endif /* ifdef(SHAKE256_H_) */