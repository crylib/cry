#include "test.h"
#include <cry/ecp.h>

static cry_ecp_grp grp;
static cry_mpi k;
static cry_ecp P, Q, Q1;

/*
 * y^2 = x^3 + 2x + 2 (mod 17)
 *
 *  Generator point = (5, 1)
 *  Generator order = 19
 */
static void simple_curve_init(cry_ecp_grp *ec)
{
    cry_mpi_init_int(&ec->a, 2);
    cry_mpi_init_int(&ec->b, 2);
    cry_mpi_init_int(&ec->p, 17);
    cry_mpi_init_int(&ec->n, 19);
    cry_mpi_init_int(&ec->g.x, 5);
    cry_mpi_init_int(&ec->g.y, 1);
    cry_mpi_init_int(&ec->g.z, 1);
}

/* Check that a point is on the curve */
static void point_check(const cry_ecp *p, const cry_ecp_grp *grp)
{
    /* Check thay y^2 = x^3 + ax + b (mod p) */
    cry_mpi v, t;

    if (cry_ecp_is_zero(p))
        return;
    cry_mpi_init_list(&v, &t, NULL);
    cry_mpi_sqr(&t, &p->x); /* x^2 */
    cry_mpi_mul(&t, &t, &p->x);         /* t = x^3 */
    cry_mpi_mul(&v, &grp->a, &p->x);    /* v = ax */
    cry_mpi_add(&v, &v, &grp->b);       /* v = ax + b */
    cry_mpi_add(&v, &v, &t);            /* v = x^3 + ax + b */
    cry_mpi_sqr(&t, &p->y);             /* t = y^2 */
    cry_mpi_sub(&v, &v, &t);
    cry_mpi_mod(&v, &v, &grp->p);
    ASSERT(cry_mpi_is_zero(&v));
    cry_mpi_clear_list(&v, &t, NULL);
}

/* Run through the whole group elements */
static void add_test(void)
{
    int i = 0;

    simple_curve_init(&grp);
    cry_ecp_init(&P);
    do {
        point_check(&P, &grp);
        cry_ecp_add(&P, &P, &grp.g, &grp);
        i++;
    } while (!cry_ecp_is_zero(&P));
    cry_ecp_grp_clear(&grp);
    cry_ecp_clear(&P);
}

static void mul_test(void)
{
    cry_mpi v;

    simple_curve_init(&grp);
    cry_ecp_init_int(&P, 9, 16);    /* p = 5g = (9, 16) */
    cry_mpi_init_int(&v, 2);

    point_check(&P, &grp);
    cry_ecp_mul(&P, &P, &v, &grp);  /* 2p = 10g = (7,11) */
    point_check(&P, &grp);

    cry_ecp_grp_clear(&grp);
    cry_ecp_clear(&P);
    cry_mpi_clear(&v);
}

static void secp192r1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_SECP192R1);

    ASSERT(res == 0);
}

static void secp224r1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_SECP224R1);

    ASSERT(res == 0);
}

static void secp256r1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_SECP256R1);

    ASSERT(res == 0);
}

static void secp384r1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_SECP384R1);

    ASSERT(res == 0);
}

static void secp521r1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_SECP521R1);

    ASSERT(res == 0);
}

static void secp192k1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_SECP192K1);

    ASSERT(res == 0);
}

static void secp224k1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_SECP224K1);

    ASSERT(res == 0);
}

static void secp256k1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_SECP256K1);

    ASSERT(res == 0);
}

static void bp256r1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_BP256R1);

    ASSERT(res == 0);
}

static void bp384r1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_BP384R1);

    ASSERT(res == 0);
}

static void bp512r1_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, CRY_ECP_GRP_BP512R1);

    ASSERT(res == 0);
}

static void bad_curve_load_test(void)
{
    int res;
    cry_ecp_grp grp;

    res = cry_ecp_grp_load(&grp, 0xFFFF);

    ASSERT(res != 0);
}

static void curve_params_load(void)
{
    RUN(secp192r1_load_test);
    RUN(secp224r1_load_test);
    RUN(secp256r1_load_test);
    RUN(secp384r1_load_test);
    RUN(secp521r1_load_test);
    RUN(secp192k1_load_test);
    RUN(secp224k1_load_test);
    RUN(secp256k1_load_test);
    RUN(bp256r1_load_test);
    RUN(bp384r1_load_test);
    RUN(bp512r1_load_test);
    RUN(bad_curve_load_test);
}


static int load_args(int argc, char *argv[])
{
    int res;

    if ((res = cry_mpi_load_str(&k, 16, argv[0])) < 0)
        return res;
    if ((res = cry_mpi_load_str(&Q.x, 16, argv[1])) < 0)
        return res;
    if ((res = cry_mpi_load_str(&Q.y, 16, argv[2])) < 0)
        return res;
    cry_mpi_set_int(&Q.z, 1);
    return res;
}

static void mul_test2(int argc, char *argv[])
{
    int res;

    ASSERT_EQ(argc, 3);
    load_args(argc, argv);

    res = cry_ecp_mul(&Q1, &P, &k, &grp);

    ASSERT_EQ(res, 0);
    ASSERT(cry_mpi_cmp(&Q.x, &Q1.x) == 0);
    ASSERT(cry_mpi_cmp(&Q.y, &Q1.y) == 0);
    ASSERT(cry_mpi_cmp(&Q.z, &Q1.z) == 0);
}

static void secp192r1_set_start_point(void)
{
    const char *x = "188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012";
    const char *y = "07192B95FFC8DA78631011ED6B24CDD573F977A11E794811";

    cry_mpi_load_str(&P.x, 16, x);
    cry_mpi_load_str(&P.y, 16, y);
    cry_mpi_set_int(&P.z, 1);
}

#ifndef SKIP_SLOW
static void secp224r1_set_start_point(void)
{
    const char *x = "B70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21";
    const char *y = "BD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34";

    cry_mpi_load_str(&P.x, 16, x);
    cry_mpi_load_str(&P.y, 16, y);
    cry_mpi_set_int(&P.z, 1);
}
#endif

#ifndef SKIP_SLOW
static void secp256r1_set_start_point(void)
{
    const char *x =
            "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296";
    const char *y =
            "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5";

    cry_mpi_load_str(&P.x, 16, x);
    cry_mpi_load_str(&P.y, 16, y);
    cry_mpi_set_int(&P.z, 1);
}
#endif

#ifndef SKIP_SLOW
static void secp384r1_set_start_point(void)
{
    const char *x =
            "AA87CA22BE8B05378EB1C71EF320AD746E1D3B628BA79B9859F741E082542A385502F25DBF55296C3A545E3872760AB7";
    const char *y =
            "3617DE4A96262C6F5D9E98BF9292DC29F8F41DBD289A147CE9DA3113B5F0B8C00A60B1CE1D7E819D7A431D7C90EA0E5F";

    cry_mpi_load_str(&P.x, 16, x);
    cry_mpi_load_str(&P.y, 16, y);
    cry_mpi_set_int(&P.z, 1);
}
#endif

#ifndef SKIP_SLOW
static void secp521r1_set_start_point(void)
{
    const char *x =
            "00C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66";
    const char *y =
            "011839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650";

    cry_mpi_load_str(&P.x, 16, x);
    cry_mpi_load_str(&P.y, 16, y);
    cry_mpi_set_int(&P.z, 1);
}
#endif

static void setup(void)
{
    cry_mpi_init(&k);
    cry_ecp_init_list(&P, &Q, &Q1, NULL);
}

static void teardown(void)
{
    cry_mpi_clear(&k);
    cry_ecp_clear_list(&P, &Q, &Q1, NULL);
}

static void dispatch(int argc, char *argv[])
{
    char *test = *argv;
    int curve = -1;

    setup();

    argv++;
    argc--;

    if (strcmp(test, "secp192r1") == 0) {
        curve = CRY_ECP_GRP_SECP192R1;
        secp192r1_set_start_point();
    } else if (strcmp(test, "secp224r1") == 0) {
#ifndef SKIP_SLOW
        curve = CRY_ECP_GRP_SECP224R1;
        secp224r1_set_start_point();
#else
        TRACE2("      (skip)\n");
#endif
    } else if (strcmp(test, "secp256r1") == 0) {
#ifndef SKIP_SLOW
        curve = CRY_ECP_GRP_SECP256R1;
        secp256r1_set_start_point();
#else
        TRACE2("      (skip)\n");
#endif
    } else if (strcmp(test, "secp384r1") == 0) {
#ifndef SKIP_SLOW
        curve = CRY_ECP_GRP_SECP384R1;
        secp384r1_set_start_point();
#else
        TRACE2("      (skip)\n");
#endif
    } else if (strcmp(test, "secp521r1") == 0) {
#ifndef SKIP_SLOW
        curve = CRY_ECP_GRP_SECP521R1;
        secp521r1_set_start_point();
#else
        TRACE2("      (skip)\n");
#endif
    } else {
        TRACE("Test '%s' not defined\n", test);
    }

    if (curve != -1) {
        cry_ecp_grp_load(&grp, curve);
        mul_test2(argc, argv);
    }

    teardown();
}

void ecp_test(void)
{
    TRACE("* ECP Algebra\n");
    curve_params_load();
    RUN(add_test);
    RUN(mul_test);
    func_test("ecp_test.data", dispatch);
    TRACE("\n");
}
