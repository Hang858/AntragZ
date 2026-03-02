#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cdt_zalcon_data.h"
#include "rng.h"

typedef __int128_t int128_t;
#define PRECISION_BITS (6 * CDT_LAYERS) 
#define PRECISION_SCALE (1ULL << PRECISION_BITS)

#define LW_LAYERS 4
#define LW_NUM_SAMPLES (1 << LW_LAYERS)

static const int64_t LW_WEIGHTS[LW_NUM_SAMPLES] = {
    10492167840LL, 10492155072LL, 10457654130LL, 10457641404LL,
    9742727280LL, 9742715424LL, 9710678835LL, 9710667018LL,
    6994778560LL, 6994770048LL, 6971769420LL, 6971760936LL,
    6495151520LL, 6495143616LL, 6473785890LL, 6473778012LL
};

static inline void FastDivMod12289(uint64_t x, uint64_t *quo, uint64_t *rem) {
    const uint64_t MAGIC = 11430799819ULL;
    const int SHIFT = 47;

#if defined(__SIZEOF_INT128__)
    unsigned __int128 prod = (unsigned __int128)x * MAGIC;
    *quo = (uint64_t)(prod >> SHIFT);
#else
    *quo = x / 12289; 
#endif
    *rem = x - (*quo * 12289);
}

#define SW_LAYERS 3
#define SW_NUM_SAMPLES 8

int32_t BaseSample(uint8_t center_idx) {
    uint64_t r = get_secure_random_u64();
    int32_t flip_mask = 0;
    
    if (center_idx > (CDT_BETA / 2)) {
        center_idx = CDT_BETA - center_idx; 
        flip_mask = 1; 
    } //利用对称性查表

    const uint64_t* cdf_row = CDT_TABLE[center_idx];
    int32_t z_index = 0;

    for (int i = 0; i < CDT_ROWS - 1; i++) {
        z_index += (r >= cdf_row[i]);
    }

    int32_t z = z_index + CDT_MIN_VAL; 
    int32_t mask = -flip_mask; 
    int32_t flipped_z = 1 - z;
    z = (flipped_z & mask) | (z & ~mask);
    return z;
}

int64_t SampleC1(int64_t c_in) {
    int64_t c0 = c_in & (CDT_BETA - 1);
    int64_t c1 = c_in >> 6;
    int32_t z = BaseSample((uint8_t)c0);
    return c1 + z;
}


int32_t SampleArbitraryCenter128(int128_t num, uint64_t den) {
    int64_t sign = 1;
    unsigned __int128 abs_num;

    if (num < 0) {
        sign = -1;
        abs_num = (unsigned __int128)(-(num + 1)) + 1; 
    } else {
        abs_num = (unsigned __int128)num;
    }

    #define PRECISION_BITS (6 * CDT_LAYERS) 
    #define PRECISION_SCALE (1ULL << PRECISION_BITS)

    uint64_t quo = (uint64_t)(abs_num / den);
    uint64_t rem = (uint64_t)(abs_num % den);

    unsigned __int128 scaled_rem = ((unsigned __int128)rem) << PRECISION_BITS;
    
    uint64_t c_frac_val = (uint64_t)(scaled_rem / den);
    uint64_t rem_frac = (uint64_t)(scaled_rem % den);


    uint64_t r_bern = get_random_range(den);
    uint64_t c_rr = c_frac_val + (r_bern < rem_frac);

    int64_t c_int = (int64_t)quo;

    c_int += (int64_t)(c_rr >> PRECISION_BITS);
    int64_t c_frac = (int64_t)(c_rr & (PRECISION_SCALE - 1));

    for (int i = 0; i < CDT_LAYERS; i++) {
        c_frac = SampleC1(c_frac);
    }
    
    int32_t result = (int32_t)(c_int + c_frac);
    return sign * result;
}

int64_t SampleLW(void) {
    int64_t z = 0;
    for (int i = 0; i < LW_NUM_SAMPLES; i++) {
        int32_t s = BaseSample(0);
        z += (int64_t)s * LW_WEIGHTS[i];
    } 
    return z;
}

int32_t SampleArbitraryCenter64(int64_t num, uint64_t den_ignored) {

    int64_t sign = 1;
    uint64_t abs_num;

    if (num < 0) {
        sign = -1;
        abs_num = (uint64_t)(-(num + 1)) + 1;
    } else {
        abs_num = (uint64_t)num;
    }

    uint64_t q1, r1;
    FastDivMod12289(abs_num, &q1, &r1);

    uint64_t integer_part = q1 >> 28;
    uint64_t frac_part_base = (q1 & 0xFFFFFFF) << 2;

    uint64_t r1_scaled = r1 * 4;
    uint64_t tail_val, tail_rem;

    FastDivMod12289(r1_scaled, &tail_val, &tail_rem);
    uint64_t c_frac_combined = frac_part_base + tail_val;
    uint64_t r_bern;
    do {
        r_bern = get_secure_random_u64() & 0x3FFF; // 0~16383
    } while (r_bern >= 12289);

    if (r_bern < tail_rem) {
        c_frac_combined += 1;
    }

    int64_t final_int = (int64_t)integer_part;
    final_int += (c_frac_combined >> 30); // PRECISION_BITS
    
    int64_t final_frac = c_frac_combined & ((1ULL << 30) - 1);

    for (int i = 0; i < CDT_LAYERS; i++) {
        final_frac = SampleC1(final_frac);
    }

    return sign * (int32_t)(final_int + final_frac);
}
