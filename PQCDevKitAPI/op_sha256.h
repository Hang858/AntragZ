#if !defined(SHA256_H_)
#define SHA256_H_
#include <stdlib.h>
#include "operator_interface.h"

#if !defined( SHA256_LEN )
#define SHA256_LEN 32    /* The length of a SHA256 hash output */
#endif

/* SHA256 context. */
typedef struct {
  unsigned long int h[8];            /* state; this is in the CPU native format */
  unsigned long Nl, Nh;              /* number of bits processed so far */
  unsigned num;                      /* number of bytes within the below */
                                     /* buffer */
  unsigned char data[64];            /* input buffer.  This is in byte vector format */
} SHA256_CTX;

static inline void SHA256_Init(SHA256_CTX *context) {
    ABORT_IF_FAIL(OP_hash_init(OP_ALG_SHA256, context, sizeof(SHA256_CTX)));
}

static inline void SHA256_Update(SHA256_CTX *context, /* context */
                  const void *input, /* input block */ 
                  unsigned int input_len)/* length of input block */
{
    ABORT_IF_FAIL(OP_hash_absorb(OP_ALG_SHA256, context, sizeof(SHA256_CTX), input, input_len));
}

static inline void SHA256_Final(unsigned char *output, SHA256_CTX *context) {
    ABORT_IF_FAIL(OP_hash_squeeze(OP_ALG_SHA256, context, sizeof(SHA256_CTX), output, SHA256_LEN));
}

static inline void SHA256_Hash(const void *input, /* input block */ 
                  unsigned int input_len, /* length of input block */
                  unsigned char *output) /* hash value */
{
    ABORT_IF_FAIL(OP_hash(OP_ALG_SHA256, 0, 32, input, input_len, 0, output));
}

#endif /* ifdef(SHA256_H_) */

