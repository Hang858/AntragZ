#ifndef PREMATRIX_H
#define PREMATRIX_H

#include <stdint.h>
#include "common.h"     
int Run_PreMatrix(
    const int8_t *f, const int8_t *g, 
    const int8_t *F, const int8_t *G, 
    PreMatrix_Output *out
);

#endif // PREMATRIX_H