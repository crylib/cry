#ifndef CRY_MPI_PVT_H_
#define CRY_MPI_PVT_H_

#include <cry/mpi.h>
#include <cry/assert.h>
#include <cry/config.h>
#include <string.h>
#include <stdint.h>
#include "../misc.h"


/* Double precision digits */

#if CRY_MPI_DIGIT_BYTES == 1
typedef uint16_t cry_mpi_dword;
#elif CRY_MPI_DIGIT_BYTES == 2
typedef uint32_t cry_mpi_dword;
#elif CRY_MPI_DIGIT_BYTES == 4
typedef uint64_t cry_mpi_dword;
#elif CRY_MPI_DIGIT_BYTES == 8
# if defined(_WIN32) || defined(__GNUC__)
typedef unsigned __int128 cry_mpi_dword;
# else
typedef uint128_t cry_mpi_dword;
# endif
#else
# error "Invalid DIGIT_MAX value"
#endif

/* Digit max value */
#ifndef CRY_MPI_QUANTUM
# include <cry/config.h>
# ifndef CRY_MPI_DEBUG_CONF
# define CRY_MPI_QUANTUM 8
# else
# define CRY_MPI_QUANTUM 1
# endif
#endif

/* Bits to digits */
#define CRY_MPI_BITS_TO_DIGS(a) \
    ((a != 0) ? (((a) - 1) / CRY_MPI_DIGIT_BITS + 1) : 0)

/* Octets to digits */
#define CRY_MPI_BYTES_TO_DIGS(a) \
    ((a != 0) ? (((a) - 1) / CRY_MPI_DIGIT_BYTES + 1) : 0)


int cry_mpi_grow(cry_mpi *a, size_t size);

int cry_mpi_shrd(cry_mpi *a, size_t n);

int cry_mpi_shld(cry_mpi *a, size_t n);

/* decrease used while the most significant digit is zero */
#define cry_mpi_adjust(a) do { \
    while ((a)->used > 0 && (a)->data[(a)->used - 1] == 0) \
        ((a)->used)--; \
    } while (0)

#define cry_mpi_set_used(a, n) do { \
    (a)->used = (n); \
    memset((a)->data, 0, (n) * CRY_MPI_DIGIT_BYTES); \
    } while (0)

#define cry_mpi_set_bit(a, bit) \
    ((a)->data[(bit) / CRY_MPI_DIGIT_BITS] |= \
        ((cry_mpi_digit)1 << ((bit) % CRY_MPI_DIGIT_BITS)))

#define cry_mpi_is_bit_set(n, bit) \
    (((n)->data[(bit) / CRY_MPI_DIGIT_BITS] & \
        ((cry_mpi_digit)1 << ((bit) % CRY_MPI_DIGIT_BITS))) != 0)


#if defined(CRY_ARCH_X86)

#define MULADD(c0, c1, c2, i, j) asm(               \
     "movl  %6, %%eax    \n\t"                      \
     "mull  %7           \n\t"                      \
     "addl  %%eax, %0    \n\t"                      \
     "adcl  %%edx, %1    \n\t"                      \
     "adcl  $0, %2       \n\t"                      \
     : "=r"(c0), "=r"(c1), "=r"(c2)                 \
     : "0"(c0), "1"(c1), "2"(c2), "m"(i), "m"(j)    \
     : "%eax", "%edx", "cc")

#elif defined(CRY_ARCH_X86_64)

#define MULADD(c0, c1, c2, i, j) asm(               \
     "movq  %6, %%rax    \n\t"                      \
     "mulq  %7           \n\t"                      \
     "addq  %%rax, %0    \n\t"                      \
     "adcq  %%rdx, %1    \n\t"                      \
     "adcq  $0, %2       \n\t"                      \
     : "=r"(c0), "=r"(c1), "=r"(c2)                 \
     : "0"(c0), "1"(c1), "2"(c2), "g"(i), "g"(j)    \
     : "%rax", "%rdx", "cc")

#elif defined(CRY_ARCH_ARM)

#define MULADD(c0, c1, c2, i, j) asm(               \
    "umull  r0, r1, %6, %7  \n\t"                   \
    "adds   %0, %0, r0      \n\t"                   \
    "adcs   %1, %1, r1      \n\t"                   \
    "adc    %2, %2, #0      \n\t"                   \
    : "=r"(c0), "=r"(c1), "=r"(c2)                  \
    : "0"(c0), "1"(c1), "2"(c2), "r"(i), "r"(j)     \
    : "r0", "r1", "cc")

#else /* ANSI C code */

#define MULADD(c0, c1, c2, i, j) do { \
   cry_mpi_dword _t; \
   _t = (cry_mpi_dword)(c0) + ((cry_mpi_dword)(i)) * ((cry_mpi_dword)(j)); \
   (c0) = (cry_mpi_digit) _t;                       \
   _t = (cry_mpi_dword)(c1) + (_t >> CRY_MPI_DIGIT_BITS); \
   (c1) = (cry_mpi_digit) _t;                       \
   (c2) += (cry_mpi_digit) (_t >> CRY_MPI_DIGIT_BITS); \
   } while (0);

#endif

#endif /* CRY_MPI_PVT_H_ */
