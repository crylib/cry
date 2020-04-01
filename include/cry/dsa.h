/**
 * @file    dsa.h
 * @brief   Digital Signature Algorithm
 */

#ifndef CRY_DSA_H_
#define CRY_DSA_H_

#include <cry/mpi.h>

struct cry_dsa_ctx {
    cry_mpi p;    /**< Prime modulus  */
    cry_mpi q;    /**< Prime modulus subgroup (prime divisor of p-1) */
    cry_mpi g;    /**< Subgroup generator */
    cry_mpi pvt;  /**< Private signing key */
    cry_mpi pub;  /**< Public verification key */
};

typedef struct cry_dsa_ctx cry_dsa_ctx;

struct cry_dsa_sig {
    cry_mpi r;  /**< Signature first component */
    cry_mpi s;  /**< Signature second component */
};

typedef struct cry_dsa_sig cry_dsa_sig;

#ifdef __cplusplus
extern "C" {
#endif

int cry_dsa_init(cry_dsa_ctx *ctx);

void cry_dsa_clear(cry_dsa_ctx *ctx);

int cry_dsa_sign(cry_dsa_ctx *ctx, cry_dsa_sig *sig,
                 const unsigned char *in, size_t len);

int cry_dsa_verify(cry_dsa_ctx *ctx, const cry_dsa_sig *sig,
                   const unsigned char *in, size_t len);

int cry_dsa_keygen(cry_dsa_ctx *ctx, unsigned int l);

#ifdef __cplusplus
}
#endif

#endif /* CRY_DSA_H_ */
