#include "mpi_pvt.h"
#include "misc.h"


int cry_mpi_mul_karatsuba(cry_mpi *r, const cry_mpi *a, const cry_mpi *b)
{
    int res;
    size_t B, hB;
    cry_mpi x0, x1, y0, y1, z0, z1, z2;

    if (cry_mpi_is_zero(a) || cry_mpi_is_zero(b)) {
        cry_mpi_zero(r);
        return 0;
    }

    /* Minimum number of digits */
    B = CRY_MIN(a->used, b->used);
    /* Divide by two */
    hB = B >> 1;
    B = hB << 1;

    /* Initialize all the temporaries */
    CRY_CHK(res = cry_mpi_init_size(&x0, hB), e0);
    CRY_CHK(res = cry_mpi_init_size(&x1, a->used - hB), e1);
    CRY_CHK(res = cry_mpi_init_size(&y0, hB), e2);
    CRY_CHK(res = cry_mpi_init_size(&y1, b->used - hB), e3);
    CRY_CHK(res = cry_mpi_init_size(&z0, B), e4);
    CRY_CHK(res = cry_mpi_init_size(&z1, B), e5);
    CRY_CHK(res = cry_mpi_init_size(&z2, B), e6);

    /* Manually set "in use" digits */
    x0.used = y0.used = hB;
    x1.used = a->used - hB;
    y1.used = b->used - hB;
    memcpy(x0.data, a->data, hB * sizeof(cry_mpi_digit));
    memcpy(y0.data, b->data, hB * sizeof(cry_mpi_digit));
    if (a->used > hB)
        memcpy(x1.data, a->data + hB, (a->used - hB) * sizeof(*x1.data));
    if (b->used > hB)
        memcpy(y1.data, b->data + hB, (b->used - hB) * sizeof(*y1.data));

    /*
     * Only need to clamp the lower words since by definition the upper words
     * x1 and y1 have a known number of digits
     */
    cry_mpi_adjust(&x0);
    cry_mpi_adjust(&y0);

    CRY_CHK(res = cry_mpi_mul_abs(&z0, &x0, &y0), e7); /* z0 = x0*y0 */
    CRY_CHK(res = cry_mpi_mul_abs(&z2, &x1, &y1), e7); /* z2 = x1*y1 */
    /* Use x0 for temporary storage */
    CRY_CHK(res = cry_mpi_add(&z1, &x1, &x0), e7);     /* z1 = x1+x0 */
    CRY_CHK(res = cry_mpi_add(&x0, &y1, &y0), e7);     /* x0 = y1+y0 */
    CRY_CHK(res = cry_mpi_mul_abs(&z1, &z1, &x0), e7); /* z1 = (x1+x0)(y1+y0) */
    CRY_CHK(res = cry_mpi_add(&x0, &z0, &z2), e7);     /* x0 = z0+z2 */
    CRY_CHK(res = cry_mpi_sub(&z1, &z1, &x0), e7); /* z1=(x1+x0)*(y1+y0)-(z0+z2) */

    CRY_CHK(res = cry_mpi_shld(&z1, hB), e7);
    CRY_CHK(res = cry_mpi_shld(&z2, B), e7);

    CRY_CHK(res = cry_mpi_add(&z1, &z0, &z1), e7);
    CRY_CHK(res = cry_mpi_add(r, &z1, &z2), e7); /* r = z2<<B + z1<<hB + z0 */

e7: cry_mpi_clear(&z2);
e6: cry_mpi_clear(&z1);
e5: cry_mpi_clear(&z0);
e4: cry_mpi_clear(&y1);
e3: cry_mpi_clear(&y0);
e2: cry_mpi_clear(&x1);
e1: cry_mpi_clear(&x0);
e0: return res;
}
