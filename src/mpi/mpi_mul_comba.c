#include "mpi_pvt.h"


int cry_mpi_mul_comba(cry_mpi *r, const cry_mpi *a, const cry_mpi *b)
{
    int res;
    size_t ix, iy, iz, tx, ty, pa;
    cry_mpi_digit c0, c1, c2, *tmpx, *tmpy;
    cry_mpi tmp, *dst;

    if (cry_mpi_is_zero(a) || cry_mpi_is_zero(b)) {
        cry_mpi_zero(r);
        return 0;
    }

    pa = a->used + b->used;
    if (r == a || r == b) {
        if ((res = cry_mpi_init_size(&tmp, pa)) != 0)
            return res;
        dst = &tmp;
    } else {
        if ((res = cry_mpi_grow(r, pa)) != 0)
            return res;
        dst = r;
    }
    cry_mpi_set_used(dst, pa);

    c1 = c2 = 0;
    for (ix = 0; ix < pa; ix++) {
        /* get offsets into the two bignums */
        ty = CRY_MIN(ix, b->used - 1);
        tx = ix - ty;

        /* setup temporary aliases */
        tmpx = a->data + tx;
        tmpy = b->data + ty;

        /*
         * this is the number of times the loop will iterate, essentially its
         * while (tx++ < a->used && ty-- >= 0)
         */
        iy = CRY_MIN(a->used - tx, ty + 1);

        /* shift accumulator right */
        c0 = c1;
        c1 = c2;
        c2 = 0;

        /* execute loop */
        for (iz = 0; iz < iy; iz++) {
            MULADD(c0, c1, c2, *tmpx, *tmpy);
            tmpx++;
            tmpy--;
        }

        /* store term */
        dst->data[ix] = c0;
    }

    if (r != dst) {
        cry_mpi_swap(r, dst);
        cry_mpi_clear(dst);
    }
    cry_mpi_adjust(r);
    return 0;
}
