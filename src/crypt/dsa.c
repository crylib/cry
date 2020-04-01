#include <cry/dsa.h>
#include "../misc.h"


#define CHK(exp) CRY_CHK(exp, e)

/*
 * Generate a non zero random secret k less than q.
 * (i.e. generate a random less than q-1 then add 1).
 *
 * k = (rand() mod (q-1)) + 1
 *
 * Is extremely important to don't reuse the same ephemeral key k.
 *
 * If the ephemeral key k is reused then DSA signature is subject
 * to a trivial attack:
 *
 * If two signatures (r1,s1) and (r2,s2) are generated using the same
 * ephemeral key k, then r = r1 = r2.
 * (An attacker can easily detect the situation)
 *
 * The private key d can be easily found:
 *
 * s1 = (M1 - d*r) * k^-1 (mod q)
 * s2 = (M2 - d*r) * k^-1 (mod q)
 * s1 - s2 = (M1 - M2) * k^-1 (mod q)
 * k = (M1 - M2) / (s1 - s2) (mod q)
 * d = (M1 - k*s1) / r (mod q)
 */
static int secret_gen(cry_mpi *k, const cry_mpi *q)
{
    int res;
    cry_mpi t, one;
    cry_mpi_digit dig = 1;

    one.sign = 0;
    one.used = 1;
    one.alloc = 0;
    one.data = &dig;

    if ((res = cry_mpi_init(&t)) != 0)
        return res;

    CHK(cry_mpi_copy(&t, q));
    CHK(cry_mpi_sub(&t, &t, &one));
    CHK(cry_mpi_rand(k, cry_mpi_count_bits(q)));
    CHK(cry_mpi_mod(k, k, &t));
    CHK(cry_mpi_add(k, k, &one));
e:  cry_mpi_clear(&t);
    return res;
}

int cry_dsa_sign(cry_dsa_ctx *ctx, cry_dsa_sig *sig,
                 const unsigned char *in, size_t len)
{
    int res;
    cry_mpi k, z;

    if ((res = cry_mpi_init_list(&k, &z, (cry_mpi *)NULL)) != 0)
        return res;

    /* z = input buffer is eventually truncated to the size of q */
    if (ctx->q.used * CRY_MPI_DIGIT_BYTES < len)
        len = ctx->q.used * CRY_MPI_DIGIT_BYTES;
    CHK(cry_mpi_load_bin(&z, in, len));

    /* k = c mod (q-1) + 1 */
    CHK(secret_gen(&k, &ctx->q));

    /* r = (g^k mod p) mod q */
    CHK(cry_mpi_mod_exp(&sig->r, &ctx->g, &k, &ctx->p));
    CHK(cry_mpi_mod(&sig->r, &sig->r, &ctx->q));

    /* s = (inv(k) * (z + d*r)) mod q */
    CHK(cry_mpi_inv(&k, &k, &ctx->q));
    CHK(cry_mpi_mul(&sig->s, &ctx->pvt, &sig->r));
    CHK(cry_mpi_add(&sig->s, &sig->s, &z));
    CHK(cry_mpi_mul(&sig->s, &sig->s, &k));
    CHK(cry_mpi_mod(&sig->s, &sig->s, &ctx->q));

e:  cry_mpi_clear_list(&k, &z, (cry_mpi *)NULL);
    return res;
}

int cry_dsa_verify(cry_dsa_ctx *ctx, const cry_dsa_sig *sig,
                   const unsigned char *in, size_t len)
{
    int res;
    cry_mpi z, w, u1, u2;

    if ((res = cry_mpi_init_list(&z, &w, &u1, &u2, (cry_mpi *)NULL)) != 0)
        return res;

    /* z = input buffer is eventually truncated to the size of q */
    if (ctx->q.used * CRY_MPI_DIGIT_BYTES < len)
        len = ctx->q.used * CRY_MPI_DIGIT_BYTES;
    CHK(cry_mpi_load_bin(&z, in, len));

    /* w = inv(s) mod q */
    CHK(cry_mpi_copy(&w, &sig->s));
    CHK(cry_mpi_inv(&w, &w, &ctx->q));

    /* u1 = (z * w) mod q */
    CHK(cry_mpi_mul(&z, &z, &w));
    CHK(cry_mpi_mod(&z, &z, &ctx->q));

    /* u2 = (r * w) mod q */
    CHK(cry_mpi_mul(&w, &sig->r, &w));
    CHK(cry_mpi_mod(&w, &w, &ctx->q));

    /* v = (((g^u1) mod p * (y^u2) mod p) mod p) mod q */
    CHK(cry_mpi_mod_exp(&u1, &ctx->g, &z, &ctx->p));
    CHK(cry_mpi_mod_exp(&u2, &ctx->pub, &w, &ctx->p));
    CHK(cry_mpi_mul(&u1, &u1, &u2));
    CHK(cry_mpi_mod(&u1, &u1, &ctx->p));
    CHK(cry_mpi_mod(&u1, &u1, &ctx->q));

    /* Check to see if v and sig match */
    res = (cry_mpi_cmp_abs(&u1, &sig->r) == 0) ? 0 : -1;
e:  cry_mpi_clear_list(&z, &w, &u1, &u2, (cry_mpi *)NULL);
    return res;
}

int cry_dsa_init(cry_dsa_ctx *ctx)
{
    return cry_mpi_init_list(&ctx->p, &ctx->q, &ctx->g, &ctx->pub, &ctx->pvt,
                             (cry_mpi *)NULL);
}

void cry_dsa_clear(cry_dsa_ctx *ctx)
{
    cry_mpi_clear_list(&ctx->p, &ctx->q, &ctx->g, &ctx->pub, &ctx->pvt,
                       (cry_mpi *)NULL);
    cry_memset(ctx, 0, sizeof(*ctx));
}

#define ITER_MAX_OUT    1024
#define ITER_MAX_IN     4096

/*
 * Assume l=8
 * 1. Generate a random prime q such that 2^150 < q < 2^160
 * 2. Generate a random integer M with 2^1023 < M < 2^1024
 * 3. Mr = M (mod 2q)
 * 4. p-1 = M-Mr (mod 2q)  (note that M-Mr is a multiple of 2q, thus is even)
 * 5. If p is prime return (p,q) else repeat from 2.
 * The entire algorithm is repeated at most ITER_MAX_OUT times.
 */
int cry_dsa_keygen(cry_dsa_ctx *ctx, unsigned int l)
{
    int res = -1;
    cry_mpi p1, q, r, one;
    unsigned int i, j;

    cry_mpi_init_list(&p1, &q, &r, &one, NULL);
    CHK(cry_mpi_set_int(&one, 1));

    for (i = 0; i < ITER_MAX_OUT; i++) {
        /* Generate a 160-bit prime q */
        CHK(cry_mpi_prime(&ctx->q, 160, NULL));
        /* Generate a prime p betweeb 512 and 1024 bits such that q | (p-1) */
        for (j = 0; j < ITER_MAX_IN; j++) {
            CHK(cry_mpi_rand(&ctx->p, 512 + l*64));
            CHK(cry_mpi_add(&r, &ctx->q, &ctx->q));
            CHK(cry_mpi_mod(&r, &ctx->p, &r));
            CHK(cry_mpi_sub(&ctx->p, &ctx->p, &r));
            CHK(cry_mpi_add(&ctx->p, &ctx->p, &one));
            if (cry_mpi_is_prime(&ctx->p))
                break;
        }
        if (j == ITER_MAX_IN)
            continue; /* Retry */
        /* Select a generator g of the unique cyclic group of order q in Zp */
        CHK(cry_mpi_sub(&p1, &ctx->p, &one));
        CHK(cry_mpi_div(&q, NULL, &p1, &ctx->q));
        CHK(cry_mpi_sub(&p1, &ctx->p, &one)); /* Required? */
        for (j = 0; j < ITER_MAX_IN; j++) {
            /* Get a random in Zp\{0} */
            CHK(cry_mpi_rand_range(&r, &p1));
            CHK(cry_mpi_add(&r, &r, &one));
            CHK(cry_mpi_mod_exp(&ctx->g, &r, &q, &ctx->p));
            if (cry_mpi_cmp_abs(&ctx->g, &one) > 0)
                break;
        }
        if (j == ITER_MAX_IN)
            continue; /* Retry */
        /* Finally compute public and public keys */
        cry_mpi_rand_range(&ctx->pvt, &q);
        cry_mpi_mod_exp(&ctx->pub, &ctx->g, &ctx->pvt, &ctx->p);
        res = 0;
        break;
    }
e:  cry_mpi_clear_list(&p1, &q, &r, &one, NULL);
    return res;
}
