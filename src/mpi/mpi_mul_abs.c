#include "mpi_pvt.h"
#include "misc.h"
#include <cry/config.h>

#define KARATSUBA_CUTOFF 64
#define TOOM3_CUTOFF     300

int cry_mpi_mul_abs(cry_mpi *r, const cry_mpi *a, const cry_mpi *b)
{
#ifdef CRY_MPI_MUL_TOOM3
    if (CRY_MIN(a->used, b->used) > TOOM3_CUTOFF)
        return cry_mpi_mul_toom3(r, a, b);
#endif

#ifdef CRY_MPI_MUL_KARATSUBA
    if (CRY_MIN(a->used, b->used) > KARATSUBA_CUTOFF)
        return cry_mpi_mul_karatsuba(r, a, b);
#endif

#ifdef CRY_MPI_MUL_COMBA
    return cry_mpi_mul_comba(r, a, b);
#else
    return cry_mpi_mul_baseline(r, a, b);
#endif
}
