#include "operator_interface.h"
#include <stdlib.h>
#include <string.h>

static uint8_t rng_state[208]; 
static int is_seeded = 0;

int init_prng(void) {
    uint8_t seed[48];
    int ret;

    ret = OP_trng(seed, sizeof(seed));
    if (ret != OP_SUCCESS) return -1;

    ret = OP_hash_init(OP_ALG_SHAKE256, rng_state, sizeof(rng_state));
    if (ret != OP_SUCCESS) return -1;


    ret = OP_hash_absorb(OP_ALG_SHAKE256, rng_state, sizeof(rng_state), seed, sizeof(seed));
    if (ret != OP_SUCCESS) return -1;

    is_seeded = 1;
    return 0;
}


void sample_fpr(double *out, size_t n) {
    if (!is_seeded) init_prng();

    size_t bytes_needed = n * 8;
    uint8_t *buf = malloc(bytes_needed); 
    
    OP_hash_squeeze(OP_ALG_SHAKE256, rng_state, sizeof(rng_state), buf, bytes_needed);

    for (size_t i = 0; i < n; i++) {
        uint64_t v;
        memcpy(&v, &buf[i * 8], 8);
        out[i] = (v >> 11) * (1.0 / 9007199254740992.0);
    }

    free(buf);
}

uint64_t get_secure_random_u64(void) {
    if (!is_seeded) init_prng();
    uint64_t v;
    OP_hash_squeeze(OP_ALG_SHAKE256, rng_state, sizeof(rng_state), (uint8_t*)&v, sizeof(v));
    return v;
}

uint64_t get_random_range(uint64_t max) {
    if (max == 0) {
        return 0; 
    }
    if (max == 1) {
        return 0; 
    }
    uint64_t limit = UINT64_MAX - (UINT64_MAX % max);

    uint64_t v;
    do {
        v = get_secure_random_u64();
    } while (v >= limit);

    return v % max;
}
