#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

// 1. 强制移除 DEBUG 定义，确保测试纯净性能
#ifdef ZITAKA_DEBUG
#undef ZITAKA_DEBUG
#endif

#include "common.h"
#include "sample.h"
#include "rng.h"

#define ITERATIONS 10  // 测试次数

// ==========================================
// CPU Cycle 计数器
// ==========================================
static inline uint64_t cpucycles(void) {
#if defined(__riscv) || defined(__riscv__)
    uint64_t cycles;
    __asm__ volatile ("rdcycle %0" : "=r" (cycles));
    return cycles;
#elif defined(__x86_64__) || defined(__i386__)
    uint32_t lo, hi;
    __asm__ volatile ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
#else
    return (uint64_t)clock();
#endif
}
int crypto_sign_keypair(PublicKey *pk, PrivateKey *sk);
int crypto_sign(uint8_t *sig, size_t *sig_len, const uint8_t *m, size_t mlen, const PrivateKey *sk);
int crypto_verify(const uint8_t *sig, size_t sig_len, const uint8_t *m, size_t mlen, const PublicKey *pk);

int main(void) {
    printf("============================================================\n");
    printf("   Zitaka / ANTRAG_Z Performance Benchmark (N=%d)   \n", ITERATIONS);
    printf("============================================================\n");

    PublicKey pk;
    PrivateKey sk;
    uint8_t msg[] = "Benchmark Message";
    size_t mlen = strlen((char*)msg);
    uint8_t sig[2048];
    size_t sig_len = sizeof(sig);
    
    uint64_t t1, t2;
    uint64_t sum_keygen = 0;
    uint64_t sum_sign_total = 0;
    uint64_t sum_offline = 0;
    uint64_t sum_verify = 0;

    int n = ANTRAG_D;
    int64_t *p1 = calloc(n, sizeof(int64_t));
    int64_t *p2 = calloc(n, sizeof(int64_t));

    printf("Warming up...\n");
    if (crypto_sign_keypair(&pk, &sk) != 1) return 1;
    
    // 预热时也要处理拒绝
    while(crypto_sign(sig, &sig_len, msg, mlen, &sk) == 0); 
    
    OfflineSamp(&sk.mat, p1, p2);
    crypto_verify(sig, sig_len, msg, mlen, &pk);

    printf("Starting Benchmark Loop...\n");

    int i;
    uint64_t total_attempts = 0;
    for (i = 0; i < ITERATIONS; i++) {
        if (i % (ITERATIONS / 10) == 0) {
            printf("Progress: %d%%\r", i * 100 / ITERATIONS);
            fflush(stdout);
        }

        // --- 1. Key Generation ---
        t1 = cpucycles();
        int kg_ret = crypto_sign_keypair(&pk, &sk);
        t2 = cpucycles();
        if (kg_ret != 1) {
            fprintf(stderr, "\nKeyGen failed (fatal) at iter %d\n", i);
            break;
        }
        sum_keygen += (t2 - t1);

        // --- 2. Offline Sample ---
        // 这是一个确定性的过程（不算拒绝），所以单独测它很稳
        t1 = cpucycles();
        OfflineSamp(&sk.mat, p1, p2);
        t2 = cpucycles();
        sum_offline += (t2 - t1);

        // --- 3. Total Signing (包含重试逻辑) ---
        t1 = cpucycles();
        int sign_ret;
        do {
            // 如果返回0(拒绝)，循环重试
            // 真实的签名时间应当包含这些重试的耗时
            total_attempts++;
            sig_len = sizeof(sig);
            sign_ret = crypto_sign(sig, &sig_len, msg, mlen, &sk);
        } while (sign_ret == 0);
        t2 = cpucycles();
        
        if (sign_ret != 1) { // 如果是 -1 (错误)，则退出
            fprintf(stderr, "\nSign error (fatal) at iter %d\n", i);
            break;
        }
        sum_sign_total += (t2 - t1);

        // --- 4. Verification ---
        t1 = cpucycles();
        int v_ret = crypto_verify(sig, sig_len, msg, mlen, &pk);
        t2 = cpucycles();
        sum_verify += (t2 - t1);
        
        if (v_ret != 1) {
            fprintf(stderr, "\nVerify failed at iter %d\n", i);
            break;
        }
    }
    
    // [修正] 防止因为意外退出导致除数错误
    // 如果 i=0 就退出了，避免除以 0
    int count = (i == 0) ? 1 : i; 

    double avg_keygen = (double)sum_keygen / count;
    double avg_sign_total = (double)sum_sign_total / count;
    double avg_offline = (double)sum_offline / count;
    double avg_verify = (double)sum_verify / count;
    
    double avg_sign_online = avg_sign_total - avg_offline;
    if (avg_sign_online < 0) avg_sign_online = 0; 

    free(p1);
    free(p2);

    printf("\n\n");
    printf("============================================================\n");
    printf("   BENCHMARK RESULTS (Average over %d runs)   \n", count);
    printf("============================================================\n");
    printf("| Operation              | Avg CPU Cycles        |\n");
    printf("|------------------------|-----------------------|\n");
    printf("| Key Generation         | %20.0f  |\n", avg_keygen);
    printf("| Verify                 | %20.0f  |\n", avg_verify);
    printf("|------------------------|-----------------------|\n");
    printf("| Sign (Total)           | %20.0f  |\n", avg_sign_total);
    printf("|   - Offline Samp       | %20.0f  |\n", avg_offline);
    printf("|   - Online Sign (Est.) | %20.0f  |\n", avg_sign_online);
    printf("============================================================\n");
    
    double avg_retries = (double)total_attempts / count;
    
    // 计算“单次尝试”的平均总耗时
    double avg_sign_per_attempt = avg_sign_total / avg_retries;

    // 真实的 Online 耗时 = 单次尝试总耗时 - 单次离线耗时
    double true_avg_online = avg_sign_per_attempt - avg_offline;
    if (true_avg_online < 0) true_avg_online = 0;

    printf("\nAnalysis (Corrected):\n");
    printf("Average Retries per Signature: %.2f\n", avg_retries);
    printf("Real Breakdown (per attempt):\n");
    printf("  - Offline: %20.0f (%.2f%%)\n", avg_offline, avg_offline/avg_sign_per_attempt*100);
    printf("  - Online : %20.0f (%.2f%%)\n", true_avg_online, true_avg_online/avg_sign_per_attempt*100);
    return 0;
}