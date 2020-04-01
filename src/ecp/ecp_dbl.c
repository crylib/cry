#include <cry/ecp.h>
#include "../misc.h"

#define CHK(exp) CRY_CHK(res = (exp), e)

int cry_ecp_dbl(cry_ecp *pr, const cry_ecp *p1, const cry_ecp_grp *grp)
{
    int res;
    cry_mpi num, den, lam;
    cry_ecp r;

    if (cry_ecp_is_zero(p1))
        return (pr != p1) ? cry_ecp_copy(pr, p1) : 0;

    if ((res = cry_mpi_init_list(&num, &den, &lam,
                                 &r.x, &r.y, &r.z, (cry_mpi *)NULL)) != 0) {
        return res;
    }

    CHK(cry_mpi_set_int(&num, 3));          /* num = 3 */
    CHK(cry_mpi_mul(&num, &num, &p1->x));   /* num = 3 * x1 */
    CHK(cry_mpi_mul(&num, &num, &p1->x));   /* num = 3 * x1^2 */
    CHK(cry_mpi_mod(&num, &num, &grp->p));
    CHK(cry_mpi_add(&num, &num, &grp->a));  /* num = 3 * x1^2 + a */
    CHK(cry_mpi_mod(&num, &num, &grp->p));

    CHK(cry_mpi_set_int(&den, 2));          /* den =  2 */
    CHK(cry_mpi_mul(&den, &den, &p1->y));   /* den =  2 * y1 */
    CHK(cry_mpi_inv(&den, &den, &grp->p));  /* den = (2 * y1)^(-1) (mod p) */
    CHK(cry_mpi_mod(&den, &den, &grp->p));

    CHK(cry_mpi_mul(&lam, &num, &den));     /* lam = num / den */
    CHK(cry_mpi_mod(&den, &den, &grp->p));

    CHK(cry_mpi_sqr(&r.x, &lam));           /* x =  lam^2 */
    CHK(cry_mpi_mod(&r.x, &r.x, &grp->p));
    CHK(cry_mpi_sub(&r.x, &r.x, &p1->x));   /* x =  lam^2 - x1 */
    CHK(cry_mpi_sub(&r.x, &r.x, &p1->x));   /* x =  lam^2 - 2*x1 */
    CHK(cry_mpi_mod(&r.x, &r.x, &grp->p));  /* x = (lam^2 - 2*x1) % p */

    CHK(cry_mpi_sub(&r.y, &p1->x, &r.x));   /* y =   x1 - x */
    CHK(cry_mpi_mul(&r.y, &r.y, &lam));     /* y =  (x1 - x) * lam */
    CHK(cry_mpi_mod(&r.y, &r.y, &grp->p));
    CHK(cry_mpi_sub(&r.y, &r.y, &p1->y));   /* y =  (x1 - x) * lam - y1 */
    CHK(cry_mpi_mod(&r.y, &r.y, &grp->p));  /* y = ((x1 - x) * lam - y1) % p */

    CHK(cry_mpi_set_int(&r.z, 1));

    cry_ecp_swap(pr, &r);
e:  cry_mpi_clear_list(&num, &den, &lam, &r.x, &r.y, &r.z, (cry_mpi *)NULL);
    return res;
}
