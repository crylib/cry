/**
 * @file    ecdsa.h
 * @brief   Elliptic Curve Digital Signature Algorithm
 */

#ifndef CRY_ECDSA_H_
#define CRY_ECDSA_H_

#include <cry/ecp.h>

struct cry_ecdsa_ctx {
    cry_ecp_grp grp;    /**< Elliptic curve group parameters */
    cry_mpi d;          /**< Private signing key */
    cry_ecp q;          /**< Public verification key */
};

typedef struct cry_ecdsa_ctx cry_ecdsa_ctx;

struct cry_ecdsa_sig {
    cry_mpi r;  /* Signature first component */
    cry_mpi s;  /* Signature second component */
};

typedef struct cry_ecdsa_sig cry_ecdsa_sig;

#ifdef __cplusplus
extern "C" {
#endif

int cry_ecdsa_sign(cry_ecdsa_ctx *ctx, cry_ecdsa_sig *sig,
                   const unsigned char *in, size_t len);

int cry_ecdsa_verify(cry_ecdsa_ctx *ctx, const cry_ecdsa_sig *sig,
                     const unsigned char *in, size_t len);

int cry_ecdsa_init(cry_ecdsa_ctx *ctx, int grp);

void cry_ecdsa_clear(cry_ecdsa_ctx *ctx);

int cry_ecdsa_keygen(cry_ecdsa_ctx *ctx);

#ifdef __cplusplus
}
#endif

#endif /* CRY_ECDSA_H_ */
