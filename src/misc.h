/*
 * A bunch of helper macros and functions for internal use.
 */

#ifndef CRY_MISC_H_
#define CRY_MISC_H_

#include <stdint.h>
#include <stdarg.h>

/** Nop to prevent some warnings */
#define NOP do {} while (0)

#define CRY_UNUSED(v) (void)(v)

/** Statically get array number of elements */
#define CRY_ARRAY_LEN(ar) (sizeof(ar)/sizeof((ar)[0]))

/** Macro used to compute the minimum of two integral values. */
#define CRY_MIN(a, b)     (((a) < (b)) ? (a) : (b))
/** Macro used to compute the maximum of two integral values. */
#define CRY_MAX(a, b)     (((a) > (b)) ? (a) : (b))

/**
 * Increments a big endian value of a give size.
 * Used to directly increment a value within a buffer.
 */
#define CRY_INCREMENT_BE(val_ptr, val_size) do { \
    size_t __i = (val_size) - 1; \
    if (++(val_ptr)[__i] == 0) \
        while (++(val_ptr)[--__i] == 0 && __i > 0); \
    } while (0)


/** Architecture independent little endian 16 bit value write. */
#define CRY_WRITE16_LE(val, dst) do { \
    ((uint8_t *) (dst))[1] = (uint8_t) (((val) >> 8) & 0xff); \
    ((uint8_t *) (dst))[0] = (uint8_t) (((val) >> 0) & 0xff); \
    } while (0)

/** Architecture independent big endian 16 bit value write. */
#define CRY_WRITE16_BE(val, dst) do { \
    ((uint8_t *) (dst))[0] = (uint8_t) (((val) >> 8) & 0xff); \
    ((uint8_t *) (dst))[1] = (uint8_t) (((val) >> 0) & 0xff); \
    } while (0)

/** Architecture independent little endian 16 bit value read. */
#define CRY_READ16_LE(val, src) \
    ((val) = ((((uint16_t) (src)[1]) << 8) | \
              (((uint16_t) (src)[0]) << 0)))

/** Architecture independent big endian 16 bit value read. */
#define CRY_READ16_BE(val, src) \
    ((val) = ((((uint16_t) (src)[0]) << 8) | \
              (((uint16_t) (src)[1]) << 0)))

/** Architecture independent little endian 32 bit value write. */
#define CRY_WRITE32_LE(val, dst) do { \
    ((uint8_t *) (dst))[3] = (uint8_t) (((val) >> 24) & 0xff); \
    ((uint8_t *) (dst))[2] = (uint8_t) (((val) >> 16) & 0xff); \
    ((uint8_t *) (dst))[1] = (uint8_t) (((val) >>  8) & 0xff); \
    ((uint8_t *) (dst))[0] = (uint8_t) (((val) >>  0) & 0xff); \
} while (0)

/** Architecture independent big endian 32 bit value write. */
#define CRY_WRITE32_BE(val, dst) do { \
    ((uint8_t *) (dst))[0] = (uint8_t) (((val) >> 24) & 0xff); \
    ((uint8_t *) (dst))[1] = (uint8_t) (((val) >> 16) & 0xff); \
    ((uint8_t *) (dst))[2] = (uint8_t) (((val) >>  8) & 0xff); \
    ((uint8_t *) (dst))[3] = (uint8_t) (((val) >>  0) & 0xff); \
} while (0)

/** Architecture independent little endian 32 bit value read. */
#define CRY_READ32_LE(val, src) \
    ((val) = ((((uint32_t) (src)[3]) << 24) | \
              (((uint32_t) (src)[2]) << 16) | \
              (((uint32_t) (src)[1]) <<  8) | \
              (((uint32_t) (src)[0]) <<  0)))

/** Architecture independent big endian 32 bit value read. */
#define CRY_READ32_BE(val, src) \
    ((val) = ((((uint32_t) (src)[0]) << 24) | \
              (((uint32_t) (src)[1]) << 16) | \
              (((uint32_t) (src)[2]) <<  8) | \
              (((uint32_t) (src)[3]) <<  0)))

/** Architecture independent little endian 64 bit value write. */
#define CRY_WRITE64_LE(val, dst) do { \
    ((uint8_t *) (dst))[7] = (uint8_t) (((val) >> 56) & 0xff); \
    ((uint8_t *) (dst))[6] = (uint8_t) (((val) >> 48) & 0xff); \
    ((uint8_t *) (dst))[5] = (uint8_t) (((val) >> 40) & 0xff); \
    ((uint8_t *) (dst))[4] = (uint8_t) (((val) >> 32) & 0xff); \
    ((uint8_t *) (dst))[3] = (uint8_t) (((val) >> 24) & 0xff); \
    ((uint8_t *) (dst))[2] = (uint8_t) (((val) >> 16) & 0xff); \
    ((uint8_t *) (dst))[1] = (uint8_t) (((val) >>  8) & 0xff); \
    ((uint8_t *) (dst))[0] = (uint8_t) (((val) >>  0) & 0xff); \
} while (0)

/** Architecture independent big endian 64 bit value write. */
#define CRY_WRITE64_BE(val, dst) do { \
    ((uint8_t *) (dst))[0] = (uint8_t) (((val) >> 56) & 0xff); \
    ((uint8_t *) (dst))[1] = (uint8_t) (((val) >> 48) & 0xff); \
    ((uint8_t *) (dst))[2] = (uint8_t) (((val) >> 40) & 0xff); \
    ((uint8_t *) (dst))[3] = (uint8_t) (((val) >> 32) & 0xff); \
    ((uint8_t *) (dst))[4] = (uint8_t) (((val) >> 24) & 0xff); \
    ((uint8_t *) (dst))[5] = (uint8_t) (((val) >> 16) & 0xff); \
    ((uint8_t *) (dst))[6] = (uint8_t) (((val) >>  8) & 0xff); \
    ((uint8_t *) (dst))[7] = (uint8_t) (((val) >>  0) & 0xff); \
} while (0)

/** Architecture independent little endian 64 bit value read. */
#define CRY_READ64_LE(val, src) \
    ((val) = ((((uint64_t) (src)[7]) << 56) | \
              (((uint64_t) (src)[6]) << 48) | \
              (((uint64_t) (src)[5]) << 40) | \
              (((uint64_t) (src)[4]) << 32) | \
              (((uint64_t) (src)[3]) << 24) | \
              (((uint64_t) (src)[2]) << 16) | \
              (((uint64_t) (src)[1]) <<  8) | \
              (((uint64_t) (src)[0]) <<  0)))

/** Architecture independent big endian 64 bit value read. */
#define CRY_READ64_BE(val, src) \
    ((val) = ((((uint64_t) (src)[0]) << 56) | \
              (((uint64_t) (src)[1]) << 48) | \
              (((uint64_t) (src)[2]) << 40) | \
              (((uint64_t) (src)[3]) << 32) | \
              (((uint64_t) (src)[4]) << 24) | \
              (((uint64_t) (src)[5]) << 16) | \
              (((uint64_t) (src)[6]) <<  8) | \
              (((uint64_t) (src)[7]) <<  0)))


/** In-place swap macro */
#define CRY_SWAP(v1, v2) do { \
    (v1) ^= (v2); \
    (v2) ^= (v1); \
    (v1) ^= (v2); \
} while (0)

/** Rotate the bits left */
#define CRY_ROTL(val, bits, size) \
    ((val) << bits | ((val) >> ((size) - (bits))))

/** Rotate the bits right */
#define CRY_ROTR(val, bits, size) \
    ((val) >> bits | ((val) << ((size) - (bits))))

#define CRY_ROTL32(val, bits) CRY_ROTL(val, bits, 32)
#define CRY_ROTR32(val, bits) CRY_ROTR(val, bits, 32)
#define CRY_ROTL64(val, bits) CRY_ROTL(val, bits, 64)
#define CRY_ROTR64(val, bits) CRY_ROTR(val, bits, 64)


/**
 * Modular inverse using Euclid algorithm.
 *
 * @param val   Value.
 * @param mod   Module.
 * @return      On success returns the positive modular inverse.
 *              If the value is not invertible returns 0.
 */
unsigned long cry_long_inv(unsigned long val, unsigned long mod);

/**
 * Secure memset.
 *
 * Implementation that should never be optimized out by the compiler.
 *
 * @param dst   Destination buffer pointer.
 * @param val   Constant filling value.
 * @param len   Destination buffer size.
 */
#define cry_memset(dst, val, len) do { \
    volatile unsigned char *__p = (unsigned char*)(dst); \
    size_t __n = (len); \
    while (__n--) \
        *__p++ = (val); \
} while (0)

/**
 * Check CRY operation result.
 * If not zero result jump to label.
 */
#define CRY_CHK(res, label) do { \
    if ((res) != 0) \
        goto label; \
} while (0)


typedef int (*cry_list_elem_init_f)(void *elem, ...);
typedef void (*cry_list_elem_clear_f)(void *elem, ...);

int cry_list_init(cry_list_elem_init_f init, cry_list_elem_clear_f clear,
                  void *first, va_list arg);

void cry_list_clear(cry_list_elem_clear_f clear,
                    void *first, va_list args);

#endif /* CRY_MISC_H_ */
