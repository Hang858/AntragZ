#ifndef SAMPLE_H
#define SAMPLE_H
#include "common.h"
void OfflineSamp(const PreMatrix_Output *key, int64_t *out_p1, int64_t *out_p2);

void OnlineSamp(const int64_t *u_hat, 
                    const int64_t *c1_num, const int64_t *c2_num,
                    int16_t *out_z1, int16_t *out_z2);

#endif