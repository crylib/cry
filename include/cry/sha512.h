/**
 * @file    sha512.h
 *
 * @brief   Secure Hash Algorithms 384 and 512 (SHA-384 and SHA-512)
 *          defined by FIPS 180-4: Secure Hash Standard (SHS).
 */

#ifndef CRY_SHA512_H_
#define CRY_SHA512_H_

#include <stdint.h>
#include <stddef.h>

#define CRY_SHA512_DIGEST_SIZE 64

struct cry_sha512_ctx {
    unsigned char data[128]; /**< Data block being processed */
    uint64_t len[2];         /**< Number of data bytes processed */
    uint64_t state[8];       /**< Intermediate digest state */
    char is384;              /**< SHA-384 */
};

typedef struct cry_sha512_ctx cry_sha512_ctx;

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Initialize SHA-512 context.
 *
 * @param ctx   SHA-512 context to be initialized.
 * @param is384 Produce a SHA-384 output.
 */
void cry_sha512_init(struct cry_sha512_ctx *ctx, char is384);

/**
 * Clear SHA-512 context.
 *
 * @param ctx   SHA-512 context to be cleared
 */
void cry_sha512_clear(cry_sha512_ctx *ctx);

/**
 * SHA-512 process buffer.
 *
 * @param ctx   SHA-512 context.
 * @param data  Buffer holding the data.
 * @param len   Length of the input data.
 */
void cry_sha512_update(struct cry_sha512_ctx *ctx, const unsigned char *data,
                       size_t len);

/**
 * Final digest.
 *
 * @param ctx   SHA-512 context
 * @param out   SHA-512 result
 */
void cry_sha512_digest(struct cry_sha512_ctx *ctx, unsigned char *out);

/**
 * SHA-512 computation of a given input.
 *
 * @param out   SHA-512 result buffer (shall be at least 64 bytes).
 * @param data  Input data buffer.
 * @param len   Input data length.
 */
void cry_sha512(unsigned char *out, const unsigned char *data, size_t len);

/**
 * SHA-384 computation of a given input.
 *
 * @param out   SHA-384 result buffer (shall be at least 48 bytes).
 * @param data  Input data buffer.
 * @param len   Input data length.
 */
void cry_sha384(unsigned char *out, const unsigned char *data, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* CRY_SHA512_H_ */
