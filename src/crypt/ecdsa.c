#include <cry/ecdsa.h>
#include "mpi/mpi_pvt.h"

#define CHK0(exp) CRY_CHK(res = (exp), e0)
#define CHK1(exp) CRY_CHK(res = (exp), e1)

/* Used by GBCS Client Test Vector */
#define NO_RAND_KEY1 \
        "\x85\x65\x82\xad\xca\x80\xfe\x30\xa3\xa6\x5d\xb0\x75\x2c\x08\xd6" \
        "\x73\xb7\x9f\x1a\x04\xfd\x8b\x64\x2a\x63\xe0\x10\xe6\x68\xb1\xef"
/* Used by GBCS Server Test Vector */
#define NO_RAND_KEY2 \
        "\x9b\xe6\xd0\x64\x88\x74\x8d\x62\x99\x5e\x0d\xed\x41\x89\x82\x56" \
        "\x62\x38\x90\x5d\xdc\xa2\xbc\x72\x1d\x44\x6e\xaa\x93\x2c\x07\x1b"
/* Non random key for tests */
#define NO_RAND_KEY NO_RAND_KEY2

int cry_ecdsa_sign(cry_ecdsa_ctx *ctx, cry_ecdsa_sig *sig,
                   const unsigned char *in, size_t len)
{
    int res;
    cry_mpi k, z;
    cry_ecp X;

    if ((res = cry_mpi_init_list(&k, &z, (cry_mpi *)NULL)) != 0)
        return res;
    CHK0(cry_ecp_init(&X));

    /* z = input buffer is eventually truncated to the size of q */
    if (ctx->grp.n.used * CRY_MPI_DIGIT_BYTES < len)
        len = ctx->grp.n.used * CRY_MPI_DIGIT_BYTES;
    CHK1(cry_mpi_load_bin(&z, in, len));

    /* This should be a random number between 0 and n-1 */
#if 0
    CHK1(cry_mpi_load_bin(&k, NO_RAND_KEY, sizeof(NO_RAND_KEY) - 1));
#else
    CHK1(cry_mpi_rand_range(&k, &ctx->grp.n));
#endif

    CHK1(cry_ecp_mul(&X, &ctx->grp.g, &k, &ctx->grp));

    /* r = x mod n */
    cry_mpi_swap(&sig->r, &X.x);

    /* s = (inv(k) * (z + r*d)) mod n */
    CHK1(cry_mpi_inv(&k, &k, &ctx->grp.n));
    CHK1(cry_mpi_mul(&sig->s, &ctx->d, &sig->r));
    CHK1(cry_mpi_add(&sig->s, &sig->s, &z));
    CHK1(cry_mpi_mul(&sig->s, &sig->s, &k));
    CHK1(cry_mpi_mod(&sig->s, &sig->s, &ctx->grp.n));

e1: cry_ecp_clear(&X);
e0: cry_mpi_clear_list(&k, &z, (cry_mpi *)NULL);
    return 0;
}

int cry_ecdsa_verify(cry_ecdsa_ctx *ctx, const cry_ecdsa_sig *sig,
                     const unsigned char *in, size_t len)
{
    int res;
    cry_mpi z, w;
    cry_ecp G, Q;

    if ((res = cry_mpi_init_list(&z, &w, (cry_mpi *)NULL)) != 0)
        return res;
    CHK0(cry_ecp_init_list(&G, &Q, (cry_ecp *)NULL));

    /* z = input buffer is eventually truncated to the size of q */
    if (ctx->grp.n.used * CRY_MPI_DIGIT_BYTES < len)
        len = ctx->grp.n.used * CRY_MPI_DIGIT_BYTES;
    CHK1(cry_mpi_load_bin(&z, in, len));

    /* w = inv(s) mod n */
    CHK1(cry_mpi_copy(&w, &sig->s));
    CHK1(cry_mpi_inv(&w, &w, &ctx->grp.n));

    /* u1 = (z * w) mod n */
    CHK1(cry_mpi_mul(&z, &z, &w));
    CHK1(cry_mpi_mod(&z, &z, &ctx->grp.n));

    /* u2 = (r * w) mod n */
    CHK1(cry_mpi_mul(&w, &sig->r, &w));
    CHK1(cry_mpi_mod(&w, &w, &ctx->grp.n));

    /* (x1, y1) = u1 * G + u2 * Q */
    CHK1(cry_ecp_mul(&G, &ctx->grp.g, &z, &ctx->grp));
    CHK1(cry_ecp_mul(&Q, &ctx->q, &w, &ctx->grp));
    CHK1(cry_ecp_add(&G, &G, &Q, &ctx->grp));

    /* r = x1 mod n */
    CHK1(cry_mpi_mod(&G.x, &G.x, &ctx->grp.n)); // inutile

    res = (cry_mpi_cmp_abs(&G.x, &sig->r) == 0) ? 0 : -1;

e1: cry_ecp_clear_list(&G, &Q, (cry_ecp *)NULL);
e0: cry_mpi_clear_list(&z, &w, &G.x, &G.y, &Q.x, &Q.y, (cry_mpi *)NULL);
    return res;
}

int cry_ecdsa_init(cry_ecdsa_ctx *ctx, int grp_id)
{
    int res;

    res = cry_mpi_init(&ctx->d);
    if (res != 0)
        return res;
    CHK0(cry_ecp_init(&ctx->q));
    if (grp_id != -1)
        CHK1(cry_ecp_grp_load(&ctx->grp, grp_id));
    else
        CHK1(cry_ecp_grp_init(&ctx->grp));
    return 0;

e1: cry_ecp_clear(&ctx->q);
e0: cry_mpi_clear(&ctx->d);
    return res;
}

void cry_ecdsa_clear(cry_ecdsa_ctx *ctx)
{
    cry_mpi_clear(&ctx->d);
    cry_ecp_clear(&ctx->q);
    cry_ecp_grp_clear(&ctx->grp);
    cry_memset(ctx, 0, sizeof(*ctx));
}

int cry_ecdsa_keygen(cry_ecdsa_ctx *ctx)
{
    int res;

    if ((res = cry_mpi_rand_range(&ctx->d, &ctx->grp.n)) == 0)
        res = cry_ecp_mul(&ctx->q, &ctx->grp.g, &ctx->d, &ctx->grp);
    return res;
}
