#if !defined(SM3_H_)
#define SM3_H_

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include "operator_interface.h"

static inline void sm3_hash(void* input, size_t input_len, uint8_t* output){
    ABORT_IF_FAIL(OP_hash(OP_ALG_SM3, 0, 32, input, input_len, 0, output));
}

typedef struct
{
    unsigned int state[8];
    unsigned int length;
    unsigned int curlen;
    unsigned char buf[64];
} SM3_STATE;

static inline void sm3_init(SM3_STATE* state){
    ABORT_IF_FAIL(OP_hash_init(OP_ALG_SM3, state, sizeof(SM3_STATE)));
}

static inline void sm3_process(SM3_STATE* state, unsigned char* buf, int len){
    ABORT_IF_FAIL(OP_hash_absorb(OP_ALG_SM3, state, sizeof(SM3_STATE), buf, len));
}

static inline void sm3_done(SM3_STATE* state, unsigned char* hash){
    ABORT_IF_FAIL(OP_hash_squeeze(OP_ALG_SM3, state, sizeof(SM3_STATE), hash, 32));
}

#endif /* ifdef(SM3_H_) */