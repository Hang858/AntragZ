#include <stdint.h>
#include <stddef.h>
#include "operator_interface.h"


int OP_vector_mul(uint16_t *z_out, const uint16_t *x_in, const uint16_t *y_in, uint16_t length, uint16_t q) {
    if (z_out == NULL || x_in == NULL || y_in == NULL) {
        return OP_FAILURE; 
    }
    if (q == 0) {
        return OP_FAILURE;
    }
    uint64_t accumulator = 0;
    for (int i = 0; i < length; ++i) {
        accumulator += (uint64_t)x_in[i] * y_in[i];
    }

    *z_out = (uint16_t)(accumulator % q);
    return OP_SUCCESS;
}
