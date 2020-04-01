#include <cry/elgamal.h>
#include "../misc.h"


#define MAX_ITER 1024

#define CHK(exp) CRY_CHK(res = exp, e)

/*
 * Generate a random secret k, such that gcd(k, p-1) = 1
 * (should be invertible modulo p-1).
 *
 * Is extremely important to don't reuse the same ephemeral key k.
 *
 * If the ephemeral key K is reused then Elgamal signature is subject
 * to a trivial attack:
 *
 * If two signatures (r1,s1) and (r2,s2) are generated using the same
 * ephemeral key k, then r = r1 = r2.
 * (An attacker can easily detect the situation)
 *
 * The private key d can be easily found:
 *
 * s1 = (M1 - d*r) * k^-1 (mod p-1)
 * s2 = (M2 - d*r) * k^-1 (mod p-1)
 * s1 - s2 = (M1 - M2) * k^-1 (mod p-1)
 * k = (M1 - M2) / (s1 - s2) (mod p-1)
 * d = (M1 - k*s1) / r (mod p-1)
 *
 * Note: if gcd(s1-s2, p-1) <> 1 then there are multiple solutions and the
 * attacker needs to check for the correct one (not a big issue for him).
 */
static int secret_gen(cry_mpi *k, const cry_mpi *t, const cry_mpi *one)
{
    int res;
    size_t iter = 0;
    cry_mpi r;

    if ((res = cry_mpi_init(&r)) != 0)
        return res;

    res = -1;
    do {
        CHK(cry_mpi_rand_range(k, t));
        CHK(cry_mpi_gcd(&r, k, t));
        if (cry_mpi_cmp(&r, one) == 0) {
            res = 0;
            break;
        }
    } while (++iter < MAX_ITER);

e:  cry_mpi_clear(&r);
    return res;
}

int cry_elgamal_sign2(cry_elgamal_ctx *ctx, cry_elgamal_sig *sign,
                      const unsigned char *in, size_t len)
{
    int res;
    cry_mpi one, k, z, t;
    cry_mpi_digit dig = 1;

    if (len > cry_mpi_count_bytes(&ctx->p))
        return -1;

    one.sign = 0;
    one.used = 1;
    one.alloc = 1;
    one.data = &dig;

    if ((res = cry_mpi_init_list(&k, &z, &t, (cry_mpi *)NULL)) != 0)
        return res;

    CHK(cry_mpi_copy(&t, &ctx->p));
    CHK(cry_mpi_sub(&t, &t, &one));

    /* gcd(k, p-1) = 1 */
    CHK(secret_gen(&k, &t, &one));

    /* r */
    CHK(cry_mpi_mod_exp(&sign->r, &ctx->g, &k, &ctx->p));

    /* s */
    CHK(cry_mpi_load_bin(&z, in, len));
    CHK(cry_mpi_mul(&sign->s, &ctx->d, &sign->r));
    CHK(cry_mpi_sub(&sign->s, &z, &sign->s));
    CHK(cry_mpi_mod(&sign->s, &sign->s, &t));
    CHK(cry_mpi_inv(&k, &k, &t));
    CHK(cry_mpi_mul(&sign->s, &sign->s, &k));
    CHK(cry_mpi_mod(&sign->s, &sign->s, &t));

e:  cry_mpi_clear_list(&k, &z, &t, (cry_mpi *)NULL);
    return res;
}

int cry_elgamal_verify2(cry_elgamal_ctx *ctx, const cry_elgamal_sig *sign,
                        const unsigned char *in, size_t len)
{
    int res;
    cry_mpi r, z;

    if (len > cry_mpi_count_bytes(&ctx->p))
        return -1;

    if ((res = cry_mpi_init_list(&r, &z, (cry_mpi *)NULL)) != 0)
        return res;

    CHK(cry_mpi_mod_exp(&r, &ctx->y, &sign->r, &ctx->p));
    CHK(cry_mpi_mod_exp(&z, &sign->r, &sign->s, &ctx->p));
    CHK(cry_mpi_mul(&r, &r, &z));
    CHK(cry_mpi_mod(&r, &r, &ctx->p));

    CHK(cry_mpi_load_bin(&z, in, len));
    CHK(cry_mpi_mod_exp(&z, &ctx->g, &z, &ctx->p));

    res = cry_mpi_cmp(&r, &z);
    if (res != 0)
        res = -1;

e:  cry_mpi_clear_list(&r, &z, (cry_mpi *)NULL);
    return res;
}

int cry_elgamal_sign(cry_elgamal_ctx *ctx, unsigned char *sign,
                     const unsigned char *in, size_t len)
{
    int res;
    cry_elgamal_sig s;
    size_t siz;

    if ((res = cry_mpi_init_list(&s.r, &s.s, (cry_mpi *)NULL)) != 0)
        return res;
    res = cry_elgamal_sign2(ctx, &s, in, len);
    if (res == 0) {
        siz = cry_mpi_count_bytes(&ctx->p);
        res = cry_mpi_store_bin(&s.r, sign, siz, 1);
        if (res == 0)
            res = cry_mpi_store_bin(&s.s, sign + siz, siz, 1);
    }
    cry_mpi_clear_list(&s.r, &s.s, (cry_mpi *)NULL);
    return res;
}

int cry_elgamal_verify(cry_elgamal_ctx *ctx, const unsigned char *sign,
                       const unsigned char *in, size_t len)
{
    int res;
    cry_elgamal_sig s;
    size_t siz;

    if ((res = cry_mpi_init_list(&s.r, &s.s, (cry_mpi *)NULL)) != 0)
        return res;
    siz = cry_mpi_count_bytes(&ctx->p);
    res = cry_mpi_load_bin(&s.r, sign, siz);
    if (res == 0) {
        res = cry_mpi_load_bin(&s.s, sign + siz, siz);
        if (res == 0)
            res = cry_elgamal_verify2(ctx, &s, in, len);
    }
    cry_mpi_clear_list(&s.r, &s.s, (cry_mpi *)NULL);
    return res;
}

int cry_elgamal_init(cry_elgamal_ctx *ctx)
{
    return cry_mpi_init_list(&ctx->p, &ctx->d, &ctx->g, &ctx->y,
                             (cry_mpi *)NULL);
}

void cry_elgamal_clear(cry_elgamal_ctx *ctx)
{
    cry_mpi_clear_list(&ctx->p, &ctx->d, &ctx->g, &ctx->y, (cry_mpi *)NULL);
    cry_memset(ctx, 0, sizeof(*ctx));
}
