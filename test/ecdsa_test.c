#include "test.h"
#include <cry/ecdsa.h>

static cry_ecdsa_ctx ecdsa;
static cry_ecdsa_sig sig;

/* key */
static const unsigned char d_dat[] = {
    0x16, 0x79, 0x7b, 0x5c, 0x0c, 0x7e, 0xd5, 0x46,
    0x1e, 0x2f, 0xf1, 0xb8, 0x8e, 0x6e, 0xaf, 0xa0,
    0x3c, 0x0f, 0x46, 0xbf, 0x07, 0x20, 0x00, 0xdf,
    0xc8, 0x30, 0xd6, 0x15,
};

static const unsigned char sha[] = {
    0x69, 0x93, 0x25, 0xd6, 0xfc, 0x8f, 0xbb, 0xb4,
    0x98, 0x1a, 0x6d, 0xed, 0x3c, 0x3a, 0x54, 0xad,
    0x2e, 0x4e, 0x3d, 0xb8, 0xa5, 0x66, 0x92, 0x01,
    0x91, 0x20, 0x64, 0xc6, 0x4e, 0x70, 0x0c, 0x13,
    0x92, 0x48, 0xcd, 0xc1,
};

static void sign_verify(void)
{
    /* Load a well known curve params */
    cry_ecp_grp_clear(&ecdsa.grp);
    ASSERT_OK(cry_ecp_grp_load(&ecdsa.grp, CRY_ECP_GRP_SECP224R1));

#if 0
    ASSERT_OK(cry_ecdsa_keygen(&ecdsa.grp, &ecdsa.d, &ecdsa.q));
#else
    cry_mpi_load_bin(&ecdsa.d, d_dat, sizeof(d_dat));
    cry_ecp_mul(&ecdsa.q, &ecdsa.grp.g, &ecdsa.d, &ecdsa.grp);
#endif

    ASSERT_OK(cry_ecdsa_sign(&ecdsa, &sig, sha, sizeof(sha)));
#if 0
    PRINT_MPI("p", &ecdsa.grp.p, 16);
    PRINT_MPI("a", &ecdsa.grp.a, 16);
    PRINT_MPI("b", &ecdsa.grp.b, 16);
    PRINT_MPI("g.x", &ecdsa.grp.g.x, 16);
    PRINT_MPI("g.y", &ecdsa.grp.g.y, 16);
    PRINT_MPI("n", &ecdsa.grp.n, 16);
    PRINT_MPI("d", &ecdsa.d, 16);
    PRINT_MPI("q.x", &ecdsa.q.x, 16);
    PRINT_MPI("q.y", &ecdsa.q.y, 16);
    PRINT_MPI("r", &sig.r, 16);
    PRINT_MPI("s", &sig.s, 16);
#endif
    ASSERT_OK(cry_ecdsa_verify(&ecdsa, &sig, sha, sizeof(sha)));
}


static void sign_verify_small(void)
{
    unsigned char v = 26;

    /* Load a well known curve params */
    cry_mpi_set_int(&ecdsa.grp.p, 17);
    cry_mpi_set_int(&ecdsa.grp.a, 2);
    cry_mpi_set_int(&ecdsa.grp.b, 10);
    cry_mpi_set_int(&ecdsa.grp.n, 19);
    cry_mpi_set_int(&ecdsa.grp.g.x, 5);
    cry_mpi_set_int(&ecdsa.grp.g.y, 1);
    cry_mpi_set_int(&ecdsa.grp.g.z, 1);

#if 0
    ASSERT_OK(cry_ecdsa_keygen(&ecdsa.grp, &ecdsa.d, &ecdsa.q));
#else
    cry_mpi_set_int(&ecdsa.d, 7);
    cry_ecp_mul(&ecdsa.q, &ecdsa.grp.g, &ecdsa.d, &ecdsa.grp);
#endif

    ASSERT_OK(cry_ecdsa_sign(&ecdsa, &sig, &v, 1));
    ASSERT_OK(cry_ecdsa_verify(&ecdsa, &sig, &v, 1));
}


static void setup(void)
{
    cry_ecdsa_init(&ecdsa, -1);
    cry_mpi_init_list(&sig.r, &sig.s, NULL);
}

static void teardown(void)
{
    cry_ecdsa_clear(&ecdsa);
    cry_mpi_clear_list(&sig.r, &sig.s, NULL);
}

#define MYRUN(name, test) run(name, test, setup, teardown)

void ecdsa_test(void)
{
    TRACE("* ECDSA\n");
    MYRUN("Sign and verify (trivially small)", sign_verify_small);
    MYRUN("Sign and verify", sign_verify);
    TRACE("\n");
}
