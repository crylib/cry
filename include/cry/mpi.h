/**
 * @file    mpi.h
 * @brief   Multiple precision integers.
 */

#ifndef CRY_MPI_H_
#define CRY_MPI_H_

#include <stddef.h>
#include <limits.h>

/* Digit max value */
#ifndef CRY_MPI_DIGIT_MAX
# include <cry/config.h>
# ifndef CRY_MPI_DEBUG_CONF
# define CRY_MPI_DIGIT_MAX ULONG_MAX
# else
# define CRY_MPI_DIGIT_MAX UCHAR_MAX
# endif
#endif

/* Digit bytes width */
#if CRY_MPI_DIGIT_MAX == 255UL
#define CRY_MPI_DIGIT_BYTES  1
#elif CRY_MPI_DIGIT_MAX == 65535UL
#define CRY_MPI_DIGIT_BYTES  2
#elif CRY_MPI_DIGIT_MAX == 4294967295UL
#define CRY_MPI_DIGIT_BYTES  4
#elif CRY_MPI_DIGIT_MAX == 18446744073709551615UL
#define CRY_MPI_DIGIT_BYTES  8
#else
# error "Invalid DIGIT_MAX value"
#endif

/** Digit type definition */
#if CRY_MPI_DIGIT_MAX == UCHAR_MAX
#define CRY_MPI_DIGIT_TYPE unsigned char
#elif CRY_MPI_DIGIT_MAX == USHRT_MAX
#define CRY_MPI_DIGIT_TYPE unsigned short
#elif CRY_MPI_DIGIT_MAX == UINT_MAX
#define CRY_MPI_DIGIT_TYPE unsigned int
#elif CRY_MPI_DIGIT_MAX == ULONG_MAX
#define CRY_MPI_DIGIT_TYPE unsigned long
#else
#error "Invalid DIGIT_MAX value"
#endif

/** Digit type definition */
typedef CRY_MPI_DIGIT_TYPE cry_mpi_digit;

/** Number of bits in one digit */
#define CRY_MPI_DIGIT_BITS  (CRY_MPI_DIGIT_BYTES << 3)

/**
 * MPI structure.
 */
struct cry_mpi {
    int sign;               /**< Non-zero if negative */
    size_t used;            /**< Number of used digits */
    size_t alloc;           /**< Allocated digits */
    cry_mpi_digit *data;    /**< Pointer to digits */
};

typedef struct cry_mpi cry_mpi;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Construction
 */

int cry_mpi_init(cry_mpi *a);

void cry_mpi_clear(cry_mpi *a);

int cry_mpi_copy(cry_mpi *d, const cry_mpi *s);

int cry_mpi_init_copy(cry_mpi *d, const cry_mpi *s);

int cry_mpi_init_size(cry_mpi *a, size_t size);

int cry_mpi_init_int(cry_mpi *a, long val);

int cry_mpi_set_int(cry_mpi *a, long val);

int cry_mpi_get_int(cry_mpi *a, long *val);

int cry_mpi_init_bin(cry_mpi *a, const void *b, size_t size);

int cry_mpi_load_bin(cry_mpi *a, const void *b, size_t size);

int cry_mpi_store_bin(const cry_mpi *a, void *b, size_t size, int pad);

int cry_mpi_init_str(cry_mpi *a, unsigned int radix, const char *s);

int cry_mpi_load_str(cry_mpi *a, unsigned int radix, const char *s);

int cry_mpi_store_str(const cry_mpi *a, unsigned int radix, char *s);

int cry_mpi_init_list(cry_mpi *a, ... /* (cry_mpi *)NULL */);

void cry_mpi_clear_list(cry_mpi *a, ... /* (cry_mpi *)NULL */);

/*
 *  Arithmetic
 */

int cry_mpi_abs(cry_mpi *d, const cry_mpi *s);

int cry_mpi_add_abs(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_sub_abs(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_mul_abs(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_div_abs(cry_mpi *q, cry_mpi *r, const cry_mpi *a,
                    const cry_mpi *b);

int cry_mpi_add(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_sub(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_mul(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_div(cry_mpi *q, cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

#define cry_mpi_mod(r, a, b) cry_mpi_div(NULL, r, a, b)

int cry_mpi_sqr(cry_mpi *r, const cry_mpi *a);

int cry_mpi_sqrt(cry_mpi *r, const cry_mpi *a);

int cry_mpi_shl(cry_mpi *c, const cry_mpi *a, size_t n);

int cry_mpi_shr(cry_mpi *c, const cry_mpi *a, size_t n);

/*
 * More advanced arithmetic algorithms
 */

int cry_mpi_exp(cry_mpi *r, const cry_mpi *b, const cry_mpi *e);

int cry_mpi_mod_exp(cry_mpi *r, const cry_mpi *b, const cry_mpi *e,
                    const cry_mpi *m);

int cry_mpi_gcd(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_lcm(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_inv(cry_mpi *r, const cry_mpi *a, const cry_mpi *m);

int cry_mpi_mul_baseline(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_mul_comba(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_mul_karatsuba(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_mul_toom3(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

int cry_mpi_sqr_baseline(cry_mpi *r, const cry_mpi *a);

int cry_mpi_sqr_comba(cry_mpi *r, const cry_mpi *a);

int cry_mpi_sqr_karatsuba(cry_mpi *r, const cry_mpi *a);

/*
 * Utilities
 */

int cry_mpi_cmp_abs(const cry_mpi *a, const cry_mpi *b);

int cry_mpi_cmp(const cry_mpi *a, const cry_mpi *b);

size_t cry_mpi_count_bits(const cry_mpi *a);

#define cry_mpi_count_bytes(a) \
    ((cry_mpi_count_bits(a) + 7) / 8)

#define cry_mpi_swap(a, b) do { \
    cry_mpi __t = *(a);         \
    *(a) = *(b);                \
    *(b) = __t;                 \
    } while (0)

#define cry_mpi_is_zero(a) \
    ((a)->used == 0)

#define cry_mpi_is_odd(a) \
    (!cry_mpi_is_zero(a) && ((a)->data[0] & 1))

#define cry_mpi_is_even(a) \
    (!cry_mpi_is_odd(a))

#define cry_mpi_is_pos(a) \
    ((a)->sign == 0)

#define cry_mpi_is_neg(a) \
    (!cry_mpi_is_pos(a))

void cry_mpi_print(const cry_mpi *a, unsigned int radix);

#define cry_mpi_zero(a) do { \
    (a)->sign = 0;           \
    (a)->used = 0;           \
} while (0)

int cry_mpi_rand(cry_mpi *a, size_t bits);

int cry_mpi_rand_range(cry_mpi *a, const cry_mpi *max);

int cry_mpi_prime(cry_mpi *a, size_t bits, unsigned int *iter);

int cry_mpi_is_prime(const cry_mpi *a);

int cry_mpi_mul_dig(cry_mpi *r, const cry_mpi *a, cry_mpi_digit b);

#ifdef __cplusplus
}
#endif

#endif /* CRY_MPI_H_ */
