/**
 * @file    arc4.h
 *
 * @brief   ARCFOUR stream cipher.
 *
 * @warning Multiple vulnerabilities have been discovered in ARC4, rendering
 *          it insecure. It is especially vulnerable when the beginning of
 *          the output keystream is not discarded, or when nonrandom or
 *          related keys are used.
 */

#ifndef CRY_ARC4_H_
#define CRY_ARC4_H_

#include <stddef.h>
#include <stdint.h>

/**
 * ARC4 cipher context.
 */
struct cry_arc4_ctx {
    unsigned char x;        /** Permutation index */
    unsigned char y;        /** Permutation index */
    unsigned char s[256];   /** Permutation table */
};

typedef struct cry_arc4_ctx cry_arc4_ctx;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Context initialization.
 *
 * @param ctx       ARC4 context.
 */
void cry_arc4_init(cry_arc4_ctx *ctx);

/**
 * Context cleanup.
 *
 * @param ctx       ARC4 context.
 */
void cry_arc4_clear(cry_arc4_ctx *ctx);

/**
 * Set the cipher key.
 *
 * @param ctx   ARC4 context.
 * @param key   Cipher key.
 * @param size  Cipher key size.
 */
void cry_arc4_key_set(cry_arc4_ctx *ctx, const unsigned char *key,
                      size_t size);

/**
 * Encryption/Decryption function.
 *
 * @param ctx   ARC4 context.
 * @param dst   Destination buffer.
 * @param src   Source buffer.
 * @param size  Number of bytes.
 */
void cry_arc4_crypt(cry_arc4_ctx *ctx, unsigned char *dst,
                    const unsigned char *src, size_t size);

#ifdef __cplusplus
}
#endif

#define cry_arc4_decrypt cry_arc4_crypt
#define cry_arc4_encrypt cry_arc4_crypt

#endif /* CRY_ARC4_H_ */
