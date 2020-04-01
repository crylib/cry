#include "mpi_pvt.h"
#include <stdlib.h>
#include <stdio.h>

int cry_mpi_grow(cry_mpi *a, size_t digs)
{
    cry_mpi_digit *tmp;

    if (digs <= a->alloc)
        return 0;

    /* ensure there are always at least MPI_PREC digits extra on top */
    digs += (CRY_MPI_QUANTUM * 2) - (digs % CRY_MPI_QUANTUM);
    /* reallocate the array */
    tmp = (cry_mpi_digit *)realloc(a->data, sizeof(cry_mpi_digit) * digs);
    if (tmp == NULL)
        return -1;
    a->data = tmp;

    /* zero excess digits */
    a->alloc = digs;

    return 0;
}

int cry_mpi_init(cry_mpi *a)
{
    /* allocate required memory and clear it */
    a->data = (cry_mpi_digit *)malloc(sizeof(cry_mpi_digit) * CRY_MPI_QUANTUM);
    if (a->data == NULL)
        return -1;

    /* finalize the initialization */
    a->used = 0;
    a->alloc = CRY_MPI_QUANTUM;
    a->sign = 0;
    return 0;
}

/*
 * Init an cry_mpi for a given size
 */
int cry_mpi_init_size(cry_mpi *a, size_t size)
{
    if (size != 0) {
        a->data = malloc(sizeof(cry_mpi_digit) * size);
        if (a->data == NULL)
            return -1;
    } else {
        a->data = NULL;
    }
    a->used = 0;
    a->alloc = size;
    a->sign = 0;
    return 0;
}

/*
 * Initialize a big number and assign it an integer value.
 */
int cry_mpi_init_int(cry_mpi *a, long val)
{
    int res;

    if ((res = cry_mpi_init(a)) == 0) {
        if ((res = cry_mpi_set_int(a, val)) < 0)
            cry_mpi_clear(a);
    }
    return res;
}

int cry_mpi_set_int(cry_mpi *a, long val)
{
    int res;
    size_t used = CRY_MPI_BYTES_TO_DIGS(sizeof(val));
    unsigned long uval;

    if ((res = cry_mpi_grow(a, used)) < 0)
        return res;

    if (val < 0) {
        a->sign = 1;
        uval = (unsigned long)-val;
    } else {
        a->sign = 0;
        uval = (unsigned long)val;
    }

    a->used = 0;
    while (uval != 0 && a->used < used) {
        a->data[a->used++] = (cry_mpi_digit)uval;
#if CRY_MPI_DIGIT_MAX != ULONG_MAX
        uval >>= CRY_MPI_DIGIT_BITS;
#endif
    }
    return res;
}

int cry_mpi_get_int(cry_mpi *a, long *val)
{
    size_t used;
    size_t i;
    long rval = 0;

    used = cry_mpi_count_bits(a);
    /* consider the sign bit */
    if (used > 8*sizeof(long) - ((a->sign == 0) ? 1 : 0))
        return -1;
    i = a->used;
    while (i-- > 0) {
#if CRY_MPI_DIGIT_MAX != ULONG_MAX
        rval <<= CRY_MPI_DIGIT_BITS;
#endif
        rval |= a->data[i];
    }
    *val = (a->sign == 0) ? rval : -rval;
    return 0;
}

void cry_mpi_clear(cry_mpi *a)
{
    if (a->data != NULL) {
        if (a->alloc != 0) /* do not free non malloced data */
            free(a->data);
        a->data = NULL;
    }
    /* reset the members to make debugging easier */
    a->alloc = 0;
    a->used = 0;
    a->sign = 0;
}

int cry_mpi_copy(cry_mpi *d, const cry_mpi *s)
{
    int res;
    unsigned int i;

    /* if dst == src do nothing */
    if (d == s)
        return 0;

    /* grow dest */
    if (d->alloc < s->used) {
        res = cry_mpi_grow(d, s->used);
        if (res != 0)
            return res;
    }

    /* zero dest and copy the parameters over */
    for (i = 0; i < s->used; i++)
        d->data[i] = s->data[i];
    d->used = s->used;
    d->sign = s->sign;
    return 0;
}

int cry_mpi_init_copy(cry_mpi *d, const cry_mpi *s)
{
    int res;

    if ((res = cry_mpi_init(d)) == 0) {
        if ((res = cry_mpi_copy(d, s)) != 0)
            cry_mpi_clear(d);
    }
    return res;
}

/*
 * Quick way to get number of bits in a word
 */
static size_t digit_bits(cry_mpi_digit d)
{
    static const unsigned char bits[256] = {
        0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    };

#if CRY_MPI_DIGIT_BYTES == 8
    if (d & 0xFFFFFFFF00000000UL) {
        if (d & 0xFFFF000000000000UL) {
            if (d & 0xFF00000000000000UL)
                return bits[d >> 56] + 56U;
            else
                return bits[d >> 48] + 48U;
        } else {
            if (d & 0x0000FF0000000000UL)
                return bits[d >> 40] + 40U;
            else
                return bits[d >> 32] + 32U;
        }
    }
#endif
#if CRY_MPI_DIGIT_BYTES >= 4
    if (d & 0xFFFF0000UL) {
        if (d & 0xFF000000UL)
            return bits[d >> 24] + 24U;
        else
            return bits[d >> 16] + 16U;
    }
#endif
#if CRY_MPI_DIGIT_BYTES >= 2
    if (d & 0xFF00UL)
        return bits[d >> 8] + 8U;
#endif
    return bits[d];
}


/*
 * Big number size in bits
 */
size_t cry_mpi_count_bits(const cry_mpi *a)
{
    size_t n = 0;

    if (a->used != 0)
        n = ((a->used-1)*CRY_MPI_DIGIT_BITS) + digit_bits(a->data[a->used-1]);
    return n;
}
