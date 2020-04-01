/**
 * @file    sha256.h
 *
 * @brief   Secure Hash Algorithms 224 and 256 (SHA-224 and SHA-256)
 *          defined by FIPS 180-4: Secure Hash Standard (SHS).
 */

#ifndef CRY_SHA256_H_
#define CRY_SHA256_H_

#include <stdint.h>
#include <stddef.h>

#define CRY_SHA256_DIGEST_SIZE 32

struct cry_sha256_ctx {
    unsigned char data[64];
    uint32_t datalen;
    uint64_t bitlen;
    uint32_t state[8];
};

typedef struct cry_sha256_ctx cry_sha256_ctx;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Initialize SHA-256 context.
 *
 * @param ctx   SHA-256 context to be initialized
 */
void cry_sha256_init(struct cry_sha256_ctx *ctx);

/**
 * Clear SHA-256 context.
 *
 * @param ctx   SHA-256 context to be cleared
 */
void cry_sha256_clear(cry_sha256_ctx *ctx);

/**
 * SHA-256 process buffer.
 *
 * @param ctx   SHA-256 context.
 * @param data  Buffer holding the data.
 * @param len   Length of the input data.
 */
void cry_sha256_update(struct cry_sha256_ctx *ctx, const unsigned char *data,
                       size_t len);

/**
 * Final digest.
 *
 * @param ctx   SHA-256 context
 * @param out   SHA-256 result
 */
void cry_sha256_digest(struct cry_sha256_ctx *ctx, unsigned char *out);

/**
 * SHA-256 computation of a given input.
 *
 * @param out   SHA-256 result buffer (shall be at least 32 bytes).
 * @param data  Input data buffer.
 * @param len   Input data length.
 */
void cry_sha256(unsigned char *out, const unsigned char *data, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* CRY_SHA256_H */
