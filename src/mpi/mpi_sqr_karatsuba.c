#include "mpi_pvt.h"

int cry_mpi_sqr_karatsuba(cry_mpi *r, const cry_mpi *a)
{
    int res;
    size_t B, hB;
    cry_mpi x0, x1, z0, z1, z2;

    if (cry_mpi_is_zero(a)) {
        cry_mpi_zero(r);
        return 0;
    }

    /* Number of digits */
    B = a->used;
    /* Divide by two */
    hB = B >> 1;
    B = hB << 1;

    /* Initialize all the temporaries */
    CRY_CHK(res = cry_mpi_init_size(&x0, hB), e0);
    CRY_CHK(res = cry_mpi_init_size(&x1, a->used - hB), e1);
    CRY_CHK(res = cry_mpi_init_size(&z0, B), e2);
    CRY_CHK(res = cry_mpi_init_size(&z1, B), e3);
    CRY_CHK(res = cry_mpi_init_size(&z2, B), e4);

    /* Manually set "in use" digits */
    x0.used = hB;
    x1.used = a->used - hB;
    memcpy(x0.data, a->data, hB * sizeof(cry_mpi_digit));
    memcpy(x1.data, a->data + hB, (a->used - hB) * sizeof(*x1.data));

    /*
     * Only need to clamp the lower word since by definition the upper word
     * x1 have a known number of digits.
     */
    cry_mpi_adjust(&x0);

    CRY_CHK(res = cry_mpi_sqr(&z0, &x0), e5);  /* z0 = x0**2 */
    CRY_CHK(res = cry_mpi_sqr(&z2, &x1), e5);  /* z2 = x1**2 */
    /* Use x0 for temporary storage */
    CRY_CHK(res = cry_mpi_add(&z1, &x1, &x0), e5);  /* z1 = x1+x0 */
    CRY_CHK(res = cry_mpi_sqr(&z1, &z1), e5);       /* z1 = (x1+x0)**2 */
    CRY_CHK(res = cry_mpi_add(&x0, &z0, &z2), e5);  /* x0 = z0+z2 */
    CRY_CHK(res = cry_mpi_sub(&z1, &z1, &x0), e5);  /* z1=(x1+x0)**2-(z0+z2) */

    CRY_CHK(res = cry_mpi_shld(&z1, hB), e5);
    CRY_CHK(res = cry_mpi_shld(&z2, B), e5);

    CRY_CHK(res = cry_mpi_add(&z1, &z0, &z1), e5);
    CRY_CHK(res = cry_mpi_add(r, &z1, &z2), e5); /* r = z2<<B + z1<<hB + z0 */

e5: cry_mpi_clear(&z2);
e4: cry_mpi_clear(&z1);
e3: cry_mpi_clear(&z0);
e2: cry_mpi_clear(&x1);
e1: cry_mpi_clear(&x0);
e0: return res;
}
