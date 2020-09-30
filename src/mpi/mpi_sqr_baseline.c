#include "mpi_pvt.h"

int cry_mpi_sqr_baseline(cry_mpi *r, const cry_mpi *a)
{
    int res;
    size_t i, j, pa;
    cry_mpi t;
    cry_mpi_dword dd, c, ch;
    cry_mpi_digit tmpx, *tmpt;

    if (cry_mpi_is_zero(a)) {
        cry_mpi_zero(r);
        return 0;
    }

    pa = a->used;
    if ((res = cry_mpi_init_size(&t, 2*pa + 1)) != 0)
        return res;

    /* default used is maximum possible size */
    cry_mpi_set_used(&t, 2*pa + 1);

    for (i = 0; i < pa; i++) {
        dd = (cry_mpi_dword)t.data[2*i] +
             (cry_mpi_dword)a->data[i] * (cry_mpi_dword)a->data[i];
        /* store lower part in the result */
        t.data[2*i] = (cry_mpi_digit)dd;
        /* get the carry */
        c = dd >> CRY_MPI_DIGIT_BITS;

        /* left hand side of a[i]*a[i] */
        tmpx = a->data[i];
        /* alias for where to store the result */
        tmpt = t.data + (2*i + 1);

        for (j = i+1; j < pa; j++) {
            /* first calculate the double product */
            dd = (cry_mpi_dword)tmpx * (cry_mpi_dword)a->data[j];
            /* check if mul by 2 will generate overflow */
            ch = (dd & (((cry_mpi_dword)1U) << (2*CRY_MPI_DIGIT_BITS - 1)))
                 >> (CRY_MPI_DIGIT_BITS - 1);
            dd <<= 1; /* mul by 2 */
            dd += ((cry_mpi_dword)*tmpt + c);
            if (ch == 0 && dd < ((cry_mpi_dword)*tmpt + c))
                ch = (cry_mpi_dword)1U << CRY_MPI_DIGIT_BITS;
            /* store lower part */
            *tmpt++ = (cry_mpi_digit)dd;
            /* get the new carry */
            c = ch | (dd >> CRY_MPI_DIGIT_BITS);
        }

        /* propagate carry */
        while (c != 0) {
            dd = (cry_mpi_dword)*tmpt + c;
            *tmpt++ = (cry_mpi_digit)dd;
            c = dd >> CRY_MPI_DIGIT_BITS;
        }
    }

    cry_mpi_adjust(&t);
    cry_mpi_swap(r, &t);
    cry_mpi_clear(&t);
    return res;
}
