#include "operator_interface.h"

typedef struct {
  unsigned long int h[8];            /* state; this is in the CPU native format */
  unsigned long Nl, Nh;              /* number of bits processed so far */
  unsigned num;                      /* number of bytes within the below */
                                     /* buffer */
  unsigned char data[64];            /* input buffer.  This is in byte vector format */
} SHA256_CTX;

void SHA256_Init (SHA256_CTX *ctx)
{
    ctx->Nl = 0;
    ctx->Nh = 0;
    ctx->num = 0;
    ctx->h[0] = 0x6A09E667UL;
    ctx->h[1] = 0xBB67AE85UL;
    ctx->h[2] = 0x3C6EF372UL;
    ctx->h[3] = 0xA54FF53AUL;
    ctx->h[4] = 0x510E527FUL;
    ctx->h[5] = 0x9B05688CUL;
    ctx->h[6] = 0x1F83D9ABUL;
    ctx->h[7] = 0x5BE0CD19UL;
}

#define SM3_IVA 0x7380166f
#define SM3_IVB 0x4914b2b9
#define SM3_IVC 0x172442d7
#define SM3_IVD 0xda8a0600
#define SM3_IVE 0xa96f30bc
#define SM3_IVF 0x163138aa
#define SM3_IVG 0xe38dee4d
#define SM3_IVH 0xb0fb0e4e

typedef struct
{
    unsigned int state[8];
    unsigned int length;
    unsigned int curlen;
    unsigned char buf[64];
} SM3_STATE;

void SM3_init(SM3_STATE *md)
{
    md->curlen = md->length = 0;
    md->state[0] = SM3_IVA;
    md->state[1] = SM3_IVB;
    md->state[2] = SM3_IVC;
    md->state[3] = SM3_IVD;
    md->state[4] = SM3_IVE;
    md->state[5] = SM3_IVF;
    md->state[6] = SM3_IVG;
    md->state[7] = SM3_IVH;
}

int OP_hash_init(uint8_t alg, void *s, int s_len) {
    switch (alg)
    {
    case OP_ALG_SHA256:
        if ((size_t)s_len < sizeof(SHA256_CTX)) {
            return OP_FAILURE;
        }
        SHA256_Init((SHA256_CTX *)s);
        return OP_SUCCESS;
    case OP_ALG_SM3:
        if ((size_t)s_len < sizeof(SM3_STATE)) {
            return OP_FAILURE;
        }
        SM3_init((SM3_STATE*)s);
        return OP_SUCCESS;
    case OP_ALG_SHAKE128:
    case OP_ALG_SHAKE256:
    case OP_ALG_SHA3_256:
    case OP_ALG_SHA3_384:
    case OP_ALG_SHA3_512:
        if (s_len < OP_STATE_SIZE_SHA3) {
            return OP_FAILURE;
        } else {
            int i;
            for(i=0;i<s_len;i+=4) {
                ((uint32_t*)s)[i]=0;
            }
        }
        break;
    default:
        return OP_FAILURE;
    }
    return 0;
}