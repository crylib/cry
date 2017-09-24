/*
 * Copyright (c) 2013-2017, Davide Galassi. All rights reserved.
 *
 * This file is part of CRY software.
 *
 * CRY is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with CRY; if not, see <http://www.gnu/licenses/>.
 */

#include "mpi_pvt.h"

/*
 * Shift left by a certain amount of digits
 */
int cry_mpi_shld(cry_mpi *a, int n)
{
    int x, res;

    if (n <= 0)
        return 0;

    /* grow to fit the new digits */
    if (a->alloc < (a->used + n)) {
        if ((res = cry_mpi_grow(a, a->used + n)) != 0)
            return res;
    }

    {
        cry_mpi_digit *top, *bottom;

        a->used += n;
        top = a->data + a->used - 1;
        bottom = (a->data + a->used - 1) - n;
        for (x = a->used - 1; x >= n; x--)
            *top-- = *bottom--;

        /* zero the lower digits */
        top = a->data;
        for (x = 0; x < n; x++)
            *top++ = 0;
    }
    return 0;
}

/*
 * Shift left by a certain bit count
 */
int cry_mpi_shl(cry_mpi *c, const cry_mpi *a, int n)
{
    cry_mpi_digit d;
    int res;

    /* copy */
    if (a != c) {
        if ((res = cry_mpi_copy(c, a)) != 0)
            return res;
    }

    if (c->alloc < (c->used + (n / CRY_MPI_DIGIT_BITS) + 1)) {
        if ((res = cry_mpi_grow(c,
                        c->used + (n / CRY_MPI_DIGIT_BITS) + 1)) != 0)
            return res;
    }

    /* shift by as many digits in the bit count */
    if (n >= CRY_MPI_DIGIT_BITS) {
        if ((res = cry_mpi_shld(c, n / CRY_MPI_DIGIT_BITS)) != 0)
            return res;
    }

    /* shift any bit count < DIGIT_BIT */
    d = n % CRY_MPI_DIGIT_BITS;
    if (d != 0) {
        cry_mpi_digit *tmpc, shift, mask, r, rr;
        size_t x;

        mask = (((cry_mpi_digit)1) << d) - 1;
        shift = CRY_MPI_DIGIT_BITS - d;
        tmpc = c->data;
        r = 0;  /* carry */
        for (x = 0; x < c->used; x++) {
            rr = (*tmpc >> shift) & mask;
            *tmpc = ((*tmpc << d) | r);
            ++tmpc;
            r = rr;
        }

        /* set final carry */
        if (r != 0)
            c->data[(c->used)++] = r;
    }
    cry_mpi_adjust(c);
    return 0;
}

