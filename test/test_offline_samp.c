/**
 * test_offline_strict.c
 * 这是一个白盒测试，用于验证 OfflineSamp 的矩阵乘法逻辑和定点数缩放逻辑是否精确。
 * 适配全静态数组定义的 PreMatrix_Output 结构体。
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include "common.h"
#include "math/poly.h"

// ==========================================
// 1. Mock (打桩) 外部依赖
// ==========================================

// [修改] 智能 Mock 控制变量
static int64_t g_mock_lw_val = 0;       // 默认返回值
static int     g_lw_call_count = 0;     // 调用计数器
static int     g_smart_mode = 0;        // 是否开启智能模式
static int64_t g_smart_val_phase1 = 0;  // 第一阶段返回值
static int64_t g_smart_val_phase2 = 0;  // 第二阶段返回值 (通常为0)
static int     g_phase1_limit = 0;      // 第一阶段调用的截止次数

// 替代 zalcon_samp.c 中的 SampleLW
int64_t SampleLW(void) {
    if (g_smart_mode) {
        int idx = g_lw_call_count++;
        if (idx < g_phase1_limit) {
            return g_smart_val_phase1;
        } else {
            return g_smart_val_phase2;
        }
    }
    return g_mock_lw_val;
}

// 替代 zalcon_samp.c 中的 SampleArbitraryCenter128
int32_t SampleArbitraryCenter128(int128_t num, uint64_t den) {
    // 确定性返回：模拟四舍五入
    int128_t half = den / 2;
    int128_t res;
    if (num >= 0)
        res = (num + half) / (int128_t)den;
    else
        res = (num - half) / (int128_t)den;
        
    return (int32_t)res;
}

// 引入待测函数声明
void OfflineSamp(const PreMatrix_Output *key, int64_t *out_p1, int64_t *out_p2);

// Mock matrix_get
static int64_t g_mock_matrix[64][ANTRAG_D];

int64_t matrix_get(const MemContext *ctx, int row, int col) {
    if (row < 64 && col < ANTRAG_D) {
        return g_mock_matrix[row][col];
    }
    return 0;
}

// ==========================================
// 2. 测试辅助工具
// ==========================================

PreMatrix_Output* alloc_empty_key(int n) {
    PreMatrix_Output *key = (PreMatrix_Output *)calloc(1, sizeof(PreMatrix_Output));
    if (!key) {
        fprintf(stderr, "Failed to allocate PreMatrix_Output\n");
        exit(1);
    }
    return key;
}

void free_key(PreMatrix_Output *key) {
    if (key) free(key);
}

// ==========================================
// 3. 测试用例
// ==========================================

void test_identity_block() {
    printf("[Test] Checking Identity Matrix Block...\n");
    int n = ANTRAG_D;
    PreMatrix_Output *key = alloc_empty_key(n);
    
    // 关闭智能模式，使用全局固定值
    g_smart_mode = 0;
    g_mock_lw_val = 1;
    
    memset(g_mock_matrix, 0, sizeof(g_mock_matrix));
    
    int64_t p1[ANTRAG_D], p2[ANTRAG_D];
    OfflineSamp(key, p1, p2);

    for(int i=0; i<n; i++) {
        assert(p1[i] == 0);
        assert(p2[i] == 0);
    }
    printf("  -> Identity block basic pass (no crash).\n");
    free_key(key);
}

void test_accumulation_logic() {
    printf("[Test] Checking Accumulation Logic (v = Ax)...\n");
    int n = ANTRAG_D;
    PreMatrix_Output *key = alloc_empty_key(n);

    // 设定：
    // x = 2^62
    // C11 = 1
    // 目标：验证 v1 = I*x + C11*x = 2*x = 2^63
    // 需要排除 MIGD 的干扰。
    
    // 分析 OfflineSamp 的 SampleLW 调用顺序：
    // 1. Id v1 (n 次)  --> 需要 x = 2^62
    // 2. Id v2 (n 次)  --> 设为 0 忽略
    // 3. Chol Col1 (n 次) --> 需要 x = 2^62 (用于 C11*x)
    // 4. Chol Col2 (n 次) --> 设为 0
    // 5. ... 后续 MIGD ... --> 全部设为 0

    g_smart_mode = 1;
    g_lw_call_count = 0;
    
    // 我们只要让前 4n 次调用覆盖 Id 和 Chol 阶段
    // 其中为了测试 v1，我们需要第 1 组(Id v1) 和 第 3 组(Chol Col1) 返回大数
    // 但简单的 smart mode 只能分两段。
    // 没关系，我们让前 3n 次都返回大数：
    // Id v1 (大), Id v2 (大), Chol Col1 (大)。
    // 这样 v1 = 2*x (符合预期), v2 = x (无所谓)
    
    g_smart_val_phase1 = 1LL << 62;
    g_smart_val_phase2 = 0;     // MIGD 阶段返回 0
    g_phase1_limit = 3 * n;     // 覆盖 Id_v1, Id_v2, Chol_Col1
    
    key->c11[0] = 1; 

    int64_t res_p1[ANTRAG_D], res_p2[ANTRAG_D];
    OfflineSamp(key, res_p1, res_p2);
    
    // v1[0] = x[0](Id) + c11[0]*x[0](Chol) = 2^62 + 2^62 = 2^63
    // p1 = 2^63 / 2^63 = 1
    printf("  -> Expected p1[0] = 1. Actual: %ld\n", res_p1[0]);
    assert(res_p1[0] == 1);
    
    printf("  -> Expected p1[1] = 1. Actual: %ld\n", res_p1[1]);
    assert(res_p1[1] == 1);

    printf("  -> Accumulation logic PASSED.\n");
    free_key(key);
}

void test_migd_b_power_logic() {
    printf("[Test] Checking MIGD b^j Logic...\n");
    int n = ANTRAG_D;
    PreMatrix_Output *key = alloc_empty_key(n);
    
    // 目标：验证 Loop j=0..K-1 中 v1 += b^j * x 的逻辑
    // 这次我们要让 Id 和 Chol 贡献为 0，只看 MIGD
    
    g_smart_mode = 1;
    g_lw_call_count = 0;
    
    // 1. Id(2n) + Chol(2n) = 4n 次调用 -> 返回 0
    // 2. MIGD j=0 (Col1 n 次) -> 返回 x
    // 3. 之后的 -> 返回 0
    
    // 为了实现“中间一段有值”，我们需要稍微 hack 一下逻辑
    // 或者利用 smart_mode 只能两段的特性，我们在 phase1 返回 0， phase2 返回 x?
    // 不行，MIGD 在后面。
    // 我们修改一下策略：让 Id 和 Chol 的 x 为 0 很简单，只要 phase1=0, limit=4n。
    // 但 phase2 会一直持续到结束。
    // MIGD j=0, j=1, j=2 都会用到 phase2 的 x。
    
    // 设 x = 2^33.
    // Id=0, Chol=0.
    // MIGD j=0: x*b^0 = x
    // MIGD j=1: x*b^1 = 32768*x
    // MIGD j=2: x*b^2 = 32768^2*x
    // Total v1 = x * (1 + b + b^2)
    // 预期结果与之前相同： p1 = 1
    
    g_smart_val_phase1 = 0;          // Id & Chol = 0
    g_smart_val_phase2 = 1LL << 33;  // MIGD = 2^33
    g_phase1_limit = 4 * n;          // 跳过前 4n 次调用
    
    int64_t p1[ANTRAG_D], p2[ANTRAG_D];
    OfflineSamp(key, p1, p2);
    
    printf("  -> With x=2^33 (MIGD only), b=2^15, k=3:\n");
    printf("     Actual p1[0] = %ld\n", p1[0]);
    
    assert(p1[0] == 1); 

    printf("  -> MIGD b^j scaling logic PASSED.\n");
    free_key(key);
}

void test_int128_overflow_safety() {
    printf("[Test] Checking INT128 Overflow Safety...\n");
    int n = ANTRAG_D;
    PreMatrix_Output *key = alloc_empty_key(n);
    
    // 复用 accumulation logic 的设置，但 x 为负
    g_smart_mode = 1;
    g_lw_call_count = 0;
    g_smart_val_phase1 = -(1LL << 62);
    g_smart_val_phase2 = 0; 
    g_phase1_limit = 3 * n;
    
    key->c11[0] = 1;
    
    int64_t p1[ANTRAG_D], p2[ANTRAG_D];
    OfflineSamp(key, p1, p2);
    
    printf("  -> Input -2^62, Accumulate twice -> -2^63\n");
    printf("  -> Expected p1[0] = -1. Actual: %ld\n", p1[0]);
    
    assert(p1[0] == -1);
    printf("  -> Negative overflow handling PASSED.\n");
    free_key(key);
}

int main() {
    printf("=== Starting Strict Tests for OfflineSamp (Smart Mock Version) ===\n");
    printf("Note: Assuming p=2^28, L=2^35, den=2^63\n");
    
    test_identity_block();
    test_accumulation_logic();
    test_migd_b_power_logic();
    test_int128_overflow_safety();
    
    printf("\n=== All Strict Tests Passed ===\n");
    return 0;
}