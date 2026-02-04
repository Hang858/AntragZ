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

#define SW_LAYERS 3
#define SW_NUM_SAMPLES 8

static const int64_t SW_WEIGHTS[SW_NUM_SAMPLES] = {
    56LL, 48LL, 42LL, 36LL, 28LL, 24LL, 21LL, 18LL
};

static const int64_t SW_COEFFS[SW_NUM_SAMPLES] = {
    1LL, -1LL, -1LL, 1LL, -1LL, 1LL, 1LL, -1LL
};

int32_t BaseSample(uint8_t center_idx) {
    uint64_t r = get_secure_random_u64();
    int32_t flip_mask = 0;
    
    if (center_idx > (CDT_BETA / 2)) {
        center_idx = CDT_BETA - center_idx; 
        flip_mask = 1; 
    }

    const uint64_t* cdf_row = CDT_TABLE[center_idx];
    int32_t z_index = 0;

    // 常数时间遍历
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

// int32_t SampleArbitraryCenter(int64_t num, uint64_t den) {

//     int64_t sign = 1;
//     uint64_t abs_num;
    
//     if (num < 0) {
//         sign = -1;
//         abs_num = (uint64_t)(-(num + 1)) + 1; 
//     } else {
//         abs_num = (uint64_t)num;
//     }

//     unsigned __int128 scaled_num = (unsigned __int128)abs_num << PRECISION_BITS;
    
//     uint64_t quo = (uint64_t)(scaled_num / den);
//     uint64_t rem = (uint64_t)(scaled_num % den);

//     uint64_t r_bern = get_random_range(den);
    
//     uint64_t c_rr = quo + (r_bern < rem);

//     int64_t c_int = (int64_t)(c_rr >> PRECISION_BITS);
//     int64_t c_frac = (int64_t)(c_rr & (PRECISION_SCALE - 1));

//     for (int i = 0; i < CDT_LAYERS; i++) {
//         c_frac = SampleC1(c_frac);
//     }
//     int32_t result = (int32_t)(c_int + c_frac);
//     return sign * result;
// }

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

int64_t SampleSW(int128_t c_num, uint64_t den) {
    int64_t z = 0;
    
    for (int i = 0; i < SW_NUM_SAMPLES; i++) {

        int128_t target_c;
        if (SW_COEFFS[i] == 1) {
            target_c = c_num;
        } else {
            target_c = -c_num;
        }
        
        int64_t s = SampleArbitraryCenter128(target_c, den);
        
        z += s * SW_WEIGHTS[i];
    }
    return z;
}