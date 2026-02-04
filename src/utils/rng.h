#ifndef RNG_H
#define RNG_H

#include <stddef.h>

void sample_fpr(double *out, size_t n);
uint64_t get_secure_random_u64(void);
uint64_t get_random_range(uint64_t max);

#endif