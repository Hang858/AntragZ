#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "operator_interface.h" 

// ============================================================
// 1. 黄金标准：朴素实现 (Copy from your original code)
// ============================================================
void poly_mul_naive_ref(int32_t *res, const int8_t *a, const int16_t *b, int n) {
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            int idx = k - i;
            int sign = 1;
            if (idx < 0) {
                idx += n;
                sign = -1; // X^N = -1
            }
            res[k] += sign * (int32_t)a[i] * b[idx];
        }
    }
}

// 声明我们要测试的目标函数 (在其他 .c 文件中定义)
// void poly_mul_int8_int16_to_32_acc_ntt(int32_t *res, const int8_t *a, const int16_t *b, int n);
// 或者直接引用你的头文件声明
extern void poly_mul_int8_int16_to_32_acc_ntt(int32_t *res, const int8_t *a, const int16_t *b, int n);

// ============================================================
// 2. 辅助工具
// ============================================================

// 生成随机 int8
int8_t rand_int8() {
    return (int8_t)(rand() % 256 - 128);
}

// 生成随机 int16
int16_t rand_int16() {
    return (int16_t)(rand() % 65536 - 32768);
}

// 生成随机 int32
int32_t rand_int32() {
    // 简单的 32 位随机数生成
    return (rand() << 16) | (rand() & 0xFFFF);
}

// ============================================================
// 3. 主测试逻辑
// ============================================================

#define TEST_ITERATIONS 100  // 测试循环次数
#define N 512

int main() {
    srand((unsigned)time(NULL)); // 随机种子

    printf("=============================================\n");
    printf("Starting Correctness Test for NTT Convolution\n");
    printf("N = %d, Iterations = %d\n", N, TEST_ITERATIONS);
    printf("=============================================\n");

    // 内存分配
    int8_t  *a = malloc(N * sizeof(int8_t));
    int16_t *b = malloc(N * sizeof(int16_t));
    
    // 两个结果数组：一个跑朴素，一个跑NTT
    int32_t *res_ref = malloc(N * sizeof(int32_t));
    int32_t *res_ntt = malloc(N * sizeof(int32_t));

    if (!a || !b || !res_ref || !res_ntt) {
        printf("Memory allocation failed!\n");
        return -1;
    }

    int total_errors = 0;

    for (int iter = 0; iter < TEST_ITERATIONS; iter++) {
        // 1. 生成随机输入数据
        for (int i = 0; i < N; i++) {
            a[i] = rand_int8();
            b[i] = rand_int16();
            
            // 为了测试 _acc (累加) 特性，我们需要给结果数组预填相同的初始值
            // 如果初始值不一致，或者累加逻辑错了，结果就会对不上
            int32_t initial_val = rand_int32() % 10000; // 随机初始噪声
            res_ref[i] = initial_val;
            res_ntt[i] = initial_val;
        }

        // 2. 运行朴素算法 (Golden Reference)
        poly_mul_naive_ref(res_ref, a, b, N);

        // 3. 运行 NTT 优化算法 (Target)
        poly_mul_int8_int16_to_32_acc_ntt(res_ntt, a, b, N);

        // 4. 结果比对
        int iter_errors = 0;
        for (int i = 0; i < N; i++) {
            if (res_ref[i] != res_ntt[i]) {
                if (iter_errors == 0) {
                    // 只打印该轮次的第一个错误，避免刷屏
                    printf("[FAIL] Iter %d, Index %d: Expected %d, Got %d (Diff: %d)\n", 
                           iter, i, res_ref[i], res_ntt[i], res_ntt[i] - res_ref[i]);
                }
                iter_errors++;
            }
        }

        if (iter_errors > 0) {
            total_errors++;
            printf("--> Iteration %d FAILED with %d mismatches.\n", iter, iter_errors);
            // 发现错误可以选择直接退出调试，或者继续跑
             break; 
        } else {
            if (iter % 10 == 0) {
                printf("[PASS] Iteration %d passed.\n", iter);
            }
        }
    }

    // 5. 总结
    printf("=============================================\n");
    if (total_errors == 0) {
        printf("✅ ALL TESTS PASSED! The NTT implementation is correct.\n");
    } else {
        printf("❌ TEST FAILED! Found errors in %d iterations.\n", total_errors);
    }
    printf("=============================================\n");

    free(a); free(b); free(res_ref); free(res_ntt);
    return total_errors;
}