#include <string.h>
#include "eigd.h"
#include "poly.h"

static const int16_t POOL_OFFSETS[24] = {
    0, 1, 2, 3, 5, 7, 9, 13, 17, 21, 29, 37, 
    45, 61, 77, 93, 125, 157, 189, 253, 317, 381, 509, 637
};

static const int16_t ROW_STRIDES[24] = {
    128, 128, 128, 64, 64, 64, 32, 32, 32, 16, 16, 16,
    8, 8, 8, 4, 4, 4, 2, 2, 2, 1, 1, 1
};

static inline int get_logic_limit(int stride) {
    return 256 / stride;
}

static inline void matrix_set(MemContext *ctx, int row, int col, int64_t val) {

    if (row >= 24) {
        if (col == 0) ctx->base_case[row - 24] = val;
        return;
    }

    if (col == 0) return;
    int stride = ROW_STRIDES[row];
    if (col % stride != 0) return;

    int logic_idx = col / stride;

    if (logic_idx >= get_logic_limit(stride)) return;
    
    if ((logic_idx & 1) == 0) return;

    int packed_idx = (logic_idx - 1) >> 1;
    int pool_idx = POOL_OFFSETS[row] + packed_idx;

    if (pool_idx < COMPRESSED_POOL_SIZE) {
        ctx->coeffs[pool_idx] = (int16_t)val;
    }
}

int64_t matrix_get(MemContext *ctx, int row, int col) {

    if (row >= 24) return (col == 0) ? ctx->base_case[row - 24] : 0;

    if (col == 0) {
        int j = row % K_VAL;
        if (j == 0) return 1;
        if (j == 1) return 32768;
        if (j == 2) return 1073741824;
        return 0;
    }

    // Case C: Coeffs
    int stride = ROW_STRIDES[row];
    if (col % stride != 0) return 0;

    int logic_idx = col / stride;
    
    if (logic_idx >= get_logic_limit(stride)) return 0;

    if ((logic_idx & 1) == 0) return 0;

    int packed_idx = (logic_idx - 1) >> 1;
    int pool_idx = POOL_OFFSETS[row] + packed_idx;

    return (int64_t)ctx->coeffs[pool_idx];
}

// 求解四平方数
static int128_t sqrt_128(int128_t n) {
    if (n < 0) return -1;
    if (n == 0) return 0;
    int128_t x = n, y = (x + 1) >> 1;
    while (y < x) { x = y; y = (x + n / x) >> 1; }
    return x;
}

static int is_square_128(int128_t n, int64_t *root) {
    if (n < 0) return 0;
    if (n == 0) { *root = 0; return 1; }
    int128_t r = sqrt_128(n);
    if (r * r == n) { *root = (int64_t)r; return 1; }
    return 0;
}

static int solve_two_squares_limited(int128_t n, int64_t *x, int64_t *y, int max_steps) {
    int128_t limit = sqrt_128(n);
    int128_t end_i = limit - max_steps;
    if (end_i < 0) end_i = 0;
    for (int128_t i = limit; i >= end_i; i--) {
        int128_t rem = n - i * i;
        int64_t r;
        if (is_square_128(rem, &r)) { *x = (int64_t)i; *y = r; return 1; }
    }
    return 0;
}

static int solve_four_squares(int128_t val, int64_t out[4]) {
    if (val < 0) return 0;
    if (val == 0) { memset(out, 0, 32); return 1; }
    int128_t limit_a = sqrt_128(val);
    for (int128_t a = limit_a; a >= 0; a--) {
        int128_t rem1 = val - a * a;
        int128_t limit_b = sqrt_128(rem1);
        int128_t b_end = limit_b - 50; 
        if (b_end < 0) b_end = 0;
        for (int128_t b = limit_b; b >= b_end; b--) {
            int128_t rem2 = rem1 - b * b;
            int64_t c, d;
            if (solve_two_squares_limited(rem2, &c, &d, 200)) {
                out[0] = (int64_t)a; out[1] = (int64_t)b; out[2] = c; out[3] = d;
                return 1;
            }
        }
    }
    return 0;
}

static int128_t calc_g_norm_sq(int64_t b, int k) {
    int128_t sum = 0, bj = 1;
    for (int i = 0; i < k; i++) { sum += bj * bj; if (i < k - 1) bj *= b; }
    return sum;
}

static void gadget_decompose(int64_t x, int64_t b, int k, int64_t *res) {
    for (int j = 0; j < k; j++) {
        int64_t r = x % b;
        if (r < 0) r += b; if (r >= b / 2) r -= b;
        res[j] = r;
        x = (x - r) / b;
    }
}



int EIGD_recursive_opt(const int64_t *input_f, int n, int current_l, int stride, 
                       int128_t d, int64_t b, int k, MemContext *ctx, int stack_offset) {
    // Base Case
    // ZITAKA_LOG("EIGD Recurse: L=%d, n=%d, d=%.1e", current_l, n, (double)d);
    if (n == 2) {
        int128_t target = d - (int128_t)input_f[0];
        if (target < 0) return 0;
        
        int64_t x[4];
        if (!solve_four_squares(target, x)) {
            ZITAKA_LOG("EIGD Fail: 4-square decomp failed for target=%.1e", (double)target);
            return 0;
        }

        ZITAKA_LOG("EIGD Base OK: %ld^2 + %ld^2 + %ld^2 + %ld^2 = %.0f", 
                   x[0], x[1], x[2], x[3], (double)target);
                   
        int base_row_start = OUTPUT_ROWS - 4;
        for (int i = 0; i < 4; i++) {
            matrix_set(ctx, base_row_start + i, 0, x[i]); 
        }
        return 1;
    }

    int half_n = n / 2;
    int128_t g_norm = calc_g_norm_sq(b, k);
    int out_row_start = (current_l - 2) * k;


    int64_t *temp_decomp = ctx->workspace;            
    int64_t *h_poly = ctx->input_stack + stack_offset;
    int64_t *poly_accum = ctx->workspace + 16;       
    int64_t *poly_tmp1  = poly_accum + n;             
    int64_t *poly_tmp2  = poly_tmp1 + n;              

    memset(poly_accum, 0, n * sizeof(int64_t));


    for (int i = 0; i < n/4; i++) { 
        int idx_low = 2 * i + 1;
        int64_t target_coeff = input_f[idx_low];
        gadget_decompose(-target_coeff, b, k, temp_decomp);
        
        for (int j = 0; j < k; j++) {
            matrix_set(ctx, out_row_start + j, idx_low * stride, temp_decomp[j]);
        }
    }
    
    for(int i=0; i<n; i++) poly_accum[i] = input_f[i];
    
    for (int j = 0; j < k; j++) {
        int row = out_row_start + j;
        for(int i=0; i<n; i++) {
            poly_tmp1[i] = matrix_get(ctx, row, i * stride);
        }
        poly_adj(poly_tmp1, poly_tmp2, n);
        poly_mul_acc_64(poly_accum, poly_tmp2, poly_tmp1, n);
    }
    
    poly_accum[0] -= (int64_t)g_norm;
    
    for (int i = 0; i < half_n; i++) {
        h_poly[i] = poly_accum[2 * i];
    }
    
    return EIGD_recursive_opt(h_poly, half_n, current_l - 1, stride * 2, d - g_norm, b, k, ctx, stack_offset + half_n);
}