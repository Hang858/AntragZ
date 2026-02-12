#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "poly_ntt.h"

#define ANTRAG_D 512
#define ANTRAG_Q 12289

// 参考实现
static inline int32_t mod_q(int32_t a) {
    int32_t r = a % ANTRAG_Q;
    return r < 0 ? r + ANTRAG_Q : r;
}

void poly_mul_mod_ref(int16_t *res, const int16_t *a, const int16_t *b) {
    memset(res, 0, ANTRAG_D * sizeof(int16_t));
    for (int i = 0; i < ANTRAG_D; i++) {
        for (int j = 0; j < ANTRAG_D; j++) {
            int k = i + j;
            int32_t val = (int32_t)a[i] * b[j];
            if (k < ANTRAG_D) res[k] = (int16_t)mod_q(res[k] + val);
            else res[k - ANTRAG_D] = (int16_t)mod_q(res[k - ANTRAG_D] - val);
        }
    }
}

int main() {
    printf("=== 512-NTT Verification (100 Tests) ===\n");
    poly_ntt_512_init();

    int16_t a[ANTRAG_D], b[ANTRAG_D];
    int16_t res_ref[ANTRAG_D], res_dut[ANTRAG_D];
    int pass_count = 0;
    int num_tests = 100;

    srand((unsigned int)time(NULL));

    for (int t = 0; t < num_tests; t++) {
        for(int i=0; i<ANTRAG_D; i++) {
            a[i] = rand() % ANTRAG_Q;
            b[i] = rand() % ANTRAG_Q;
        }

        poly_mul_mod_ref(res_ref, a, b);
        poly_ntt_512_mul(res_dut, a, b);

        int match = 1;
        for (int i = 0; i < ANTRAG_D; i++) {
            if (res_ref[i] != res_dut[i]) {
                match = 0;
                printf("[FAIL] Test %d at index %d: Ref=%d, Dut=%d\n", t, i, res_ref[i], res_dut[i]);
                break;
            }
        }
        if (match) pass_count++;
    }

    printf("\nPassed %d/%d tests.\n", pass_count, num_tests);
    return (pass_count == num_tests) ? 0 : 1;
}