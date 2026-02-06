#ifndef BENCH_H
#define BENCH_H

#include <stdint.h>
#include <stdio.h>

// 获取 CPU 周期数 (x86_64)
static inline uint64_t cpucycles(void) {
    uint64_t lo, hi;
    __asm__ volatile ("rdtsc" : "=a" (lo), "=d" (hi));
    return (hi << 32) | lo;
}

// 简单的性能统计结构
typedef struct {
    uint64_t start;
    uint64_t accumulated;
    uint64_t count;
} bench_timer;

#ifdef ZITAKA_DEBUG
    #define BENCH_INIT(name) bench_timer t_##name = {0}
    #define BENCH_START(name) t_##name.start = cpucycles()
    #define BENCH_STOP(name) t_##name.accumulated += (cpucycles() - t_##name.start); t_##name.count++
    #define BENCH_PRINT(name, label) \
        printf("[BENCH] %-20s: %lu cycles (avg)\n", label, t_##name.count > 0 ? t_##name.accumulated / t_##name.count : 0)
#else
    #define BENCH_INIT(name)
    #define BENCH_START(name)
    #define BENCH_STOP(name)
    #define BENCH_PRINT(name, label)
#endif

#endif