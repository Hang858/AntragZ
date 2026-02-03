#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../src/migd.h" 


static const int16_t POOL_OFFSETS_TEST[24] = {
    0, 1, 2, 3, 5, 7, 9, 13, 17, 21, 29, 37, 
    45, 61, 77, 93, 125, 157, 189, 253, 317, 381, 509, 637
};

static const int16_t ROW_STRIDES_TEST[24] = {
    128, 128, 128, 64, 64, 64, 32, 32, 32, 16, 16, 16,
    8, 8, 8, 4, 4, 4, 2, 2, 2, 1, 1, 1
};

// 辅助：获取逻辑索引上限
static inline int get_logic_limit_test(int stride) {
    return 256 / stride;
}

// 解压单行数据
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
    int limit = get_logic_limit_test(stride);

    for (int logic_idx = 1; logic_idx < limit; logic_idx += 2) {
        int packed_idx = (logic_idx - 1) >> 1;
        poly_out[logic_idx * stride] = ctx->coeffs[pool_base + packed_idx];
    }

    // 加上偏移量 b^j
    int j = row % K_VAL;
    int64_t b_pow = 0;
    if (j == 0) b_pow = 1;
    else if (j == 1) b_pow = 32768;
    else if (j == 2) b_pow = 1073741824;
    
    poly_out[0] += b_pow;
}

// ==========================================
// 2. 多项式运算 (验证专用)
// ==========================================
typedef __int128_t int128_t;

void poly_adj(const int64_t *src, int64_t *dst) {
    dst[0] = src[0];
    for (int i = 1; i < N_MAX; i++) dst[i] = -src[N_MAX - i];
}

// 累加乘积到 int128 数组: target += a * b
void poly_mul_acc_128(int128_t *target, const int64_t *a, const int64_t *b) {
    for (int i = 0; i < N_MAX; i++) {
        for (int j = 0; j < N_MAX; j++) {
            int k = i + j;
            int128_t val = (int128_t)a[i] * (int128_t)b[j];
            if (k < N_MAX) target[k] += val;
            else           target[k - N_MAX] -= val;
        }
    }
}

// 计算 ||g||^2
int128_t get_g_norm_sq() {
    int128_t sum = 0, b = 1;
    for(int i=0; i<3; i++) { sum += b*b; if(i<2) b *= 32768; }
    return sum;
}

// ==========================================
// 3. 核心验证函数
// ==========================================

// 验证 EIGD 分解: sum(row * row*) == target_d - f_poly
int verify_eigd_component(const char* name, const MemContext *ctx, const int64_t *f_poly, int128_t target_d) {
    printf("  [Check] Verifying %s component...\n", name);
    
    int128_t *sum = (int128_t *)calloc(N_MAX, sizeof(int128_t));
    int64_t *row = (int64_t *)malloc(N_MAX * sizeof(int64_t));
    int64_t *adj = (int64_t *)malloc(N_MAX * sizeof(int64_t));

    // 1. 计算向量模长平方
    for (int r = 0; r < OUTPUT_ROWS; r++) {
        expand_row_poly(ctx, r, row);
        poly_adj(row, adj);
        poly_mul_acc_128(sum, row, adj);
    }

    // 2. 验证等式: sum[i] + f_poly[i] == (i==0 ? target_d : 0)
    int errs = 0;
    
    // 常数项
    int128_t check0 = sum[0] + f_poly[0];
    if (check0 != target_d) {
        printf("    [FAIL] %s: Const term mismatch! Diff: %lld\n", name, (long long)(check0 - target_d));
        errs++;
    }

    // 高次项
    for (int i = 1; i < N_MAX; i++) {
        int128_t val = sum[i] + f_poly[i];
        if (val != 0) {
            if (errs < 5) printf("    [FAIL] %s: Index %d mismatch! Val: %lld\n", name, i, (long long)val);
            errs++;
        }
    }

    free(sum); free(row); free(adj);
    return (errs == 0);
}

// 验证 Gadget 分解: sum(b^j * c_j) == -sigma_12
int verify_gadget_component(const int64_t *c_coeffs, const int64_t *sigma_12) {
    printf("  [Check] Verifying c component (Gadget)...\n");
    
    int128_t *sum = (int128_t *)calloc(N_MAX, sizeof(int128_t));
    int64_t b_pow = 1;

    for (int j = 0; j < K_VAL; j++) {
        const int64_t *c_j = c_coeffs + (j * N_MAX);
        for(int i=0; i<N_MAX; i++) {
            sum[i] += (int128_t)c_j[i] * b_pow;
        }
        b_pow *= 32768;
    }

    int errs = 0;
    for (int i = 0; i < N_MAX; i++) {
        // Target is -sigma_12
        int128_t target = -sigma_12[i];
        if (sum[i] != target) {
            if (errs < 5) printf("    [FAIL] Gadget Index %d mismatch! Calc: %lld, Target: %lld\n", i, (long long)sum[i], (long long)target);
            errs++;
        }
    }
    
    free(sum);
    return (errs == 0);
}

// ==========================================
// 4. Main
// ==========================================
void gen_self_adjoint(int64_t *f) {
    memset(f, 0, N_MAX * sizeof(int64_t));
    f[0] = (rand() % 1000) - 500;
    for (int i = 1; i < N_MAX/2; i++) {
        int64_t val = (rand() % 500) - 250;
        f[i] = val; f[N_MAX - i] = -val;
    }
}

void gen_random(int64_t *f) {
    for(int i=0; i<N_MAX; i++) f[i] = (rand() % 500) - 250;
}

int main() {
    srand(time(NULL));
    printf("=== Strict MIGD Verification Test (Optimized Storage) ===\n");

    // 1. 准备输入
    int64_t *s11 = malloc(N_MAX * 8);
    int64_t *s12 = malloc(N_MAX * 8);
    int64_t *s22 = malloc(N_MAX * 8);
    
    gen_self_adjoint(s11);
    gen_random(s12);
    gen_self_adjoint(s22);

    int128_t d = get_g_norm_sq() * 9 * 2; // 估算一个足够大的 d
    int128_t d_prime = d - get_g_norm_sq();

    // 2. 运行 MIGD
    MIGD_Output out;
    memset(&out, 0, sizeof(MIGD_Output)); // 重要：清零

    printf("[Step 1] Running MIGD...\n");
    clock_t start = clock();
    int success = Run_MIGD(s11, s12, s22, d, 32768, &out);
    clock_t end = clock();

    if (!success) {
        printf("[FATAL] MIGD Failed.\n");
        return 1;
    }
    printf("  -> Finished in %.2f ms.\n", (double)(end - start) * 1000.0 / CLOCKS_PER_SEC);

    // 3. 验证步骤
    printf("\n[Step 2] Verifying Components...\n");
    
    int ok_x = verify_eigd_component("x (Sigma_11)", &out.x_ctx, s11, d_prime);
    int ok_c = verify_gadget_component(out.c_coeffs, s12);
    
    // 为了验证 y，我们需要手动计算 Schur 补
    int64_t *Pi = calloc(N_MAX, 8);
    int64_t *adj = malloc(N_MAX * 8);
    // Init Pi = Sigma_22
    for(int i=0; i<N_MAX; i++) Pi[i] = s22[i];
    // Add sum(c* c)
    for(int j=0; j<K_VAL; j++) {
        int64_t *c = out.c_coeffs + j*N_MAX;
        poly_adj(c, adj);
        // 这里只是为了验证，用个临时的简单乘法即可
        // 注意：这里需要累加到 Pi (int64)。如果在验证代码里溢出，说明 Sigma 选得太大。
        // 测试时输入比较小，应该没事。
        int128_t *tmp_mul = calloc(N_MAX, 16);
        poly_mul_acc_128(tmp_mul, adj, c);
        for(int k=0; k<N_MAX; k++) Pi[k] += (int64_t)tmp_mul[k];
        free(tmp_mul);
    }
    free(adj);

    int ok_y = verify_eigd_component("y (Schur)", &out.y_ctx, Pi, d_prime);
    
    free(Pi); free(s11); free(s12); free(s22);

    if (ok_x && ok_c && ok_y) {
        printf("\n=== [PASS] All MIGD properties verified! ===\n");
        printf("AA* = dI - Sigma holds.\n");
        return 0;
    } else {
        printf("\n=== [FAIL] Verification failed! ===\n");
        return 1;
    }
}