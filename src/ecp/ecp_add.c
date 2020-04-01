#include <cry/ecp.h>
#include "../misc.h"

#define CHK0(exp) CRY_CHK(res = (exp), e0)
#define CHK1(exp) CRY_CHK(res = (exp), e1)

int cry_ecp_add(cry_ecp *pr, const cry_ecp *p1, const cry_ecp *p2,
                const cry_ecp_grp *grp)
{
    int res;
    cry_mpi num, den, lam;
    cry_ecp r;

    /* Check if one of the two points is the point at infinity. */
    if (cry_ecp_is_zero(p1))
        return (pr != p2) ? cry_ecp_copy(pr, p2) : 0;
    if (cry_ecp_is_zero(p2))
        return (pr != p1) ? cry_ecp_copy(pr, p1) : 0;

    if ((res = cry_mpi_init_list(&num, &den, (cry_mpi *)NULL)) != 0)
        return res;

    CHK0(cry_mpi_sub(&num, &p2->y, &p1->y)); /* num = y2 - y1 */
    CHK0(cry_mpi_mod(&num, &num, &grp->p));
    CHK0(cry_mpi_sub(&den, &p2->x, &p1->x)); /* den = x2 - x1 */
    CHK0(cry_mpi_mod(&den, &den, &grp->p));

    /*
     * If x1=x2 then:
     * - if y1=y2 obviously p1=p2
     * - else p1 is p2 inverse
     * In the second case we assume that both points are in the same EC.
     * Because the EC graph is symmetric wrt the x-axis, by construction
     * x1=x2 and y1<>y2 implies that y1=-y2, thus p1 is the p2 inverse.
     */
    if (cry_mpi_is_zero(&den)) {
        if (cry_mpi_is_zero(&num))
            res = cry_ecp_dbl(pr, p1, grp);
        else
            cry_ecp_set_zero(pr);
        goto e0;
    }
    /* Continue with lambda computation */
    CHK0(cry_mpi_init_list(&lam, &r.x, &r.y, &r.z, (cry_mpi *)NULL));
    CHK1(cry_mpi_inv(&den, &den, &grp->p));  /* den^(-1) (mod p) */
    CHK1(cry_mpi_mul(&lam, &num, &den));     /* lam = num / den */
    CHK1(cry_mpi_mod(&lam, &lam, &grp->p));

    CHK1(cry_mpi_sqr(&r.x, &lam));           /* x =  lam^2 */
    CHK1(cry_mpi_mod(&r.x, &r.x, &grp->p));
    CHK1(cry_mpi_sub(&r.x, &r.x, &p1->x));   /* x =  lam^2 - x1 */
    CHK1(cry_mpi_sub(&r.x, &r.x, &p2->x));   /* x =  lam^2 - x1 - x2 */
    CHK1(cry_mpi_mod(&r.x, &r.x, &grp->p));  /* x = (lam^2 - x1 - x2) % p */

    CHK1(cry_mpi_sub(&r.y, &p1->x, &r.x));   /* y =   x1 - x */
    CHK1(cry_mpi_mul(&r.y, &r.y, &lam));     /* y =  (x1 - x) * lam */
    CHK1(cry_mpi_mod(&r.y, &r.y, &grp->p));
    CHK1(cry_mpi_sub(&r.y, &r.y, &p1->y));   /* y =  (x1 - x) * lam - y1 */
    CHK1(cry_mpi_mod(&r.y, &r.y, &grp->p));  /* y = ((x1 - x) * lam - y1) % p */

    CHK1(cry_mpi_set_int(&r.z, 1));

    cry_ecp_swap(pr, &r);
e1: cry_mpi_clear_list(&lam, &r.x, &r.y, &r.z, (cry_mpi *)NULL);
e0: cry_mpi_clear_list(&num, &den, (cry_mpi *)NULL);
    return res;
}
