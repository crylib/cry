/**
 * @file    elgamal.h
 * @brief   Elgamal digital signature.
 */

#ifndef CRY_ELGAMAL_H_
#define CRY_ELGAMAL_H_

#include <cry/mpi.h>

/**
 * Elgamal context
 */
struct cry_elgamal_ctx {
    cry_mpi p;  /**< Prime */
    cry_mpi g;  /**< Group generator */
    cry_mpi y;  /**< g^x mod p */
    cry_mpi d;  /**< Secret exponent */
};

typedef struct cry_elgamal_ctx cry_elgamal_ctx;

struct cry_elgamal_sig {
    cry_mpi r;
    cry_mpi s;
};

typedef struct cry_elgamal_sig cry_elgamal_sig;

#ifdef __cplusplus
extern "C" {
#endif


int cry_elgamal_init(cry_elgamal_ctx *ctx);

void cry_elgamal_clear(cry_elgamal_ctx *ctx);

int cry_elgamal_sign(cry_elgamal_ctx *ctx, unsigned char *sign,
                     const unsigned char *in, size_t len);

int cry_elgamal_verify(cry_elgamal_ctx *ctx, const unsigned char *sign,
                       const unsigned char *in, size_t len);

int cry_elgamal_sign2(cry_elgamal_ctx *ctx, cry_elgamal_sig *sign,
                      const unsigned char *in, size_t len);

int cry_elgamal_verify2(cry_elgamal_ctx *ctx, const cry_elgamal_sig *sign,
                        const unsigned char *in, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* CRY_ELGAMAL_H_ */
