#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "eigd.h"

// 复用压缩表进行解压测试
static const int16_t POOL_OFFSETS_TEST[24] = {
    0, 1, 2, 3, 5, 7, 9, 13, 17, 21, 29, 37, 
    45, 61, 77, 93, 125, 157, 189, 253, 317, 381, 509, 637
};

static const int16_t ROW_STRIDES_TEST[24] = {
    128, 128, 128, 64, 64, 64, 32, 32, 32, 16, 16, 16,
    8, 8, 8, 4, 4, 4, 2, 2, 2, 1, 1, 1
};

// ==========================================
// 解压逻辑 (Reconstruction)
// ==========================================
void expand_row_poly(const MemContext *ctx, int row, int64_t *poly_out) {
    memset(poly_out, 0, N_MAX * sizeof(int64_t));

    // Base Case (int64)
    if (row >= 24) {
        poly_out[0] = ctx->base_case[row - 24];
        return;
    }

    // Recursive Rows (int16)
    int stride = ROW_STRIDES_TEST[row];
    int pool_base = POOL_OFFSETS_TEST[row];
    int limit = N_MAX / 2; // 只读取前一半逻辑空间

    for (int logic_idx = 1; logic_idx * stride < limit; logic_idx += 2) {
        int packed_idx = (logic_idx - 1) >> 1;
        poly_out[logic_idx * stride] = ctx->coeffs[pool_base + packed_idx];
    }

    // 加上偏移量
    int j = row % K_VAL;
    int64_t b_pow = 0;
    if (j == 0) b_pow = 1;
    else if (j == 1) b_pow = 32768;
    else if (j == 2) b_pow = 1073741824;
    
    poly_out[0] += b_pow;
}

// ==========================================
// 验证逻辑
// ==========================================
void poly_adj(const int64_t *src, int64_t *dst) {
    dst[0] = src[0];
    for (int i = 1; i < N_MAX; i++) dst[i] = -src[N_MAX - i];
}

void poly_mul_accum(int128_t *target, const int64_t *a, const int64_t *b) {
    for (int i = 0; i < N_MAX; i++) {
        for (int j = 0; j < N_MAX; j++) {
            int k = i + j;
            int128_t val = (int128_t)a[i] * (int128_t)b[j];
            if (k < N_MAX) target[k] += val; else target[k - N_MAX] -= val;
        }
    }
}

int verify_compressed(MemContext *ctx, int64_t *f_input, int128_t d) {
    printf("[Verify] Reconstructing matrix and checking AA* = d - f...\n");
    int128_t *sum_poly = (int128_t *)calloc(N_MAX, sizeof(int128_t));
    int64_t *row_poly = (int64_t *)malloc(N_MAX * sizeof(int64_t));
    int64_t *row_adj = (int64_t *)malloc(N_MAX * sizeof(int64_t));

    for (int r = 0; r < OUTPUT_ROWS; r++) {
        expand_row_poly(ctx, r, row_poly);
        poly_adj(row_poly, row_adj);
        poly_mul_accum(sum_poly, row_poly, row_adj);
    }

    int errors = 0;
    int128_t check_0 = sum_poly[0] + f_input[0];
    if (check_0 != d) {
        printf("  [FAIL] Index 0 mismatch! Diff: %lld\n", (long long)(check_0 - d));
        errors++;
    }

    for (int i = 1; i < N_MAX; i++) {
        int128_t check_i = sum_poly[i] + f_input[i];
        if (check_i != 0) {
            if (errors < 5) printf("  [FAIL] Index %d mismatch! Val: %lld\n", i, (long long)check_i);
            errors++;
        }
    }

    free(sum_poly); free(row_poly); free(row_adj);
    return (errors == 0);
}

// 辅助函数
void gen_self_adjoint(int64_t *f) {
    memset(f, 0, N_MAX * sizeof(int64_t));
    f[0] = (rand() % 1000) - 500; 
    for (int i = 1; i < N_MAX / 2; i++) {
        int64_t val = (rand() % 1000) - 500;
        f[i] = val; f[N_MAX - i] = -val; 
    }
}

int128_t calc_d(int64_t b, int k, int l) {
    int128_t sum = 0, bj = 1;
    for(int i=0; i<k; i++) { sum += bj*bj; if(i<k-1) bj*=b; }
    int128_t total = sum * l;
    return total + (sum >> 1);
}

int main() {
    srand(time(NULL));
    printf("=== EIGD Optimized Storage Test ===\n");

    MemContext ctx;
    memset(&ctx, 0, sizeof(MemContext)); // 清零很重要！
    
    int64_t *f = (int64_t *)malloc(N_MAX * sizeof(int64_t));
    gen_self_adjoint(f);
    int128_t d = calc_d(32768, K_VAL, L_VAL);

    clock_t start = clock();
    int success = EIGD_recursive_opt(f, N_MAX, L_VAL, 1, d, 32768, K_VAL, &ctx, 0);
    clock_t end = clock();

    if (success) {
        printf("[SUCCESS] EIGD Finished in %.2f ms.\n", (double)(end - start) * 1000.0 / CLOCKS_PER_SEC);
        printf("Memory Usage (Key Data):\n");
        printf("  Coeffs: %lu bytes\n", sizeof(ctx.coeffs));
        printf("  Base:   %lu bytes\n", sizeof(ctx.base_case));
        printf("  Total:  < 2.1 KB\n");
        
        if (verify_compressed(&ctx, f, d)) {
            printf("[PASS] Verification Successful!\n");
        } else {
            printf("[FAIL] Verification Failed.\n");
        }
    } else {
        printf("[FAIL] Algorithm returned 0.\n");
    }

    free(f);
    return 0;
}