#include <cry/ecp.h>

int cry_ecp_init(cry_ecp *p)
{
    int res;

    res = cry_mpi_init_list(&p->x, &p->y, &p->z, (cry_mpi *)NULL);
    if (res == 0)
        cry_ecp_set_zero(p);
    return res;
}

void cry_ecp_clear(cry_ecp *p)
{
    cry_mpi_clear_list(&p->x, &p->y, &p->z, (cry_mpi *)NULL);
}

int cry_ecp_init_int(cry_ecp *p, long x, long y)
{
    int res;

    if ((res = cry_ecp_init(p)) == 0) {
        if ((res = cry_mpi_set_int(&p->x, x)) == 0) {
            if ((res = cry_mpi_set_int(&p->y, y)) == 0)
                res = cry_mpi_set_int(&p->z, 1);
        }
        if (res != 0)
            cry_ecp_clear(p);
    }
    return res;
}

int cry_ecp_copy(cry_ecp *d, const cry_ecp *s)
{
    int res;

    if ((res = cry_mpi_copy(&d->x, &s->x)) == 0) {
        if ((res = cry_mpi_copy(&d->y, &s->y)) == 0)
            res = cry_mpi_copy(&d->z, &s->z);
    }
    return res;
}


int cry_ecp_grp_init(cry_ecp_grp *grp)
{
    int res;

    res = cry_mpi_init_list(&grp->p, &grp->a, &grp->b, &grp->n,
                            &grp->g.x, &grp->g.y, &grp->g.z,
                            (cry_mpi *)NULL);
    if (res == 0)
        cry_ecp_set_zero(&grp->g);
    return res;
}

void cry_ecp_grp_clear(cry_ecp_grp *grp)
{
    cry_mpi_clear_list(&grp->p, &grp->a, &grp->b, &grp->n,
                       &grp->g.x, &grp->g.y, &grp->g.z,
                       (cry_mpi *)NULL);
}
