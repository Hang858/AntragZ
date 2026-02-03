#ifndef MIGD_H
#define MIGD_H

#include "common.h"


int Run_MIGD(
    const int64_t *sigma_11, 
    const int64_t *sigma_12, 
    const int64_t *sigma_22,
    int128_t d,
    int64_t b, 
    MIGD_Output *out_A
);

#endif // MIGD_H