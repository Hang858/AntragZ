#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "hash.h"
#include <stdlib.h>
#include "encode.h"
#include "poly.h"
#include "bench.h"
static inline int16_t center_mod(int32_t a) {
    int32_t q = ANTRAG_Q;
    int32_t r = a % q;
    if (r > q / 2) r -= q;
    if (r < -q / 2) r += q;
    return (int16_t)r;
}

int crypto_verify(const uint8_t *sig, size_t sig_len, const uint8_t *m, size_t mlen, const PublicKey *pk) {
    int n = ANTRAG_D;
    BENCH_INIT(total);
    BENCH_INIT(hash);
    BENCH_INIT(decompress);
    BENCH_INIT(poly_mul);
    
    BENCH_START(total);
    
    int16_t s1[ANTRAG_D];
    int16_t c_hash[ANTRAG_D];
    int16_t tmp_s2[ANTRAG_D];

    if (sig_len < 41) return 0;
    BENCH_START(hash);
    hash_to_point(c_hash, sig, 40, m, mlen);
    BENCH_STOP(hash);
    BENCH_START(decompress);
    size_t bytes_read = decompress_sig(s1, sig + 40, sig_len - 40);
    BENCH_STOP(decompress);
    if (bytes_read == 0) return 0;
    BENCH_START(poly_mul);
    poly_mul_mod(tmp_s2, s1, pk->h);
    // poly_mul_mod_opt(tmp_s2, s1, pk->h);
    BENCH_STOP(poly_mul);
    int64_t norm_sq = 0;
    int64_t beta_sq = 34715664LL;

    for (int i = 0; i < n; i++) {
        int16_t s2_i = center_mod((int32_t)tmp_s2[i] + c_hash[i]);
        
        norm_sq += (int64_t)s1[i] * s1[i];
        norm_sq += (int64_t)s2_i * s2_i;
    }
    BENCH_STOP(total);
    printf("=== Verify Performance (Cycles) ===\n");
    BENCH_PRINT(hash,       "Hashing");
    BENCH_PRINT(decompress, "Decompression");
    BENCH_PRINT(poly_mul,   "Poly Mul (Naive)");
    BENCH_PRINT(total,      "Total Verify");
    printf("===================================\n");
    if (norm_sq <= beta_sq) {
        ZITAKA_LOG("VERIFY OK. Norm: %ld <= Limit: %ld", norm_sq, beta_sq);
        return 1;
    } else {
        ZITAKA_LOG("VERIFY FAIL. Norm: %ld > Limit: %ld (Diff: %ld)", 
                   norm_sq, beta_sq, norm_sq - beta_sq);
        return 0;
    }
}