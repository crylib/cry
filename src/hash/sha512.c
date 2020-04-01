#include <cry/sha512.h>
#include <string.h>
#include "../misc.h"


#define S0(x) (CRY_ROTR64(x, 1) ^ CRY_ROTR64(x, 8) ^  ((x) >> 7))
#define S1(x) (CRY_ROTR64(x,19) ^ CRY_ROTR64(x,61) ^  ((x) >> 6))

#define S2(x) (CRY_ROTR64(x,28) ^ CRY_ROTR64(x,34) ^ CRY_ROTR64(x,39))
#define S3(x) (CRY_ROTR64(x,14) ^ CRY_ROTR64(x,18) ^ CRY_ROTR64(x,41))

#define F0(x,y,z) (((x) & (y)) | ((z) & ((x) | (y))))
#define F1(x,y,z) ((z) ^ ((x) & ((y) ^ (z))))

#define P(a, b, c, d, e, f, g, h, x, K) do { \
    t1 = (h) + S3(e) + F1((e),(f),(g)) + (K) + (x); \
    t2 = S2(a) + F0((a),(b),(c));   \
    (d) += t1; \
    (h) = t1 + t2; \
} while (0)


/*
 * Round constants
 */
static const uint64_t K[80] = {
    0x428A2F98D728AE22, 0x7137449123EF65CD,
    0xB5C0FBCFEC4D3B2F, 0xE9B5DBA58189DBBC,
    0x3956C25BF348B538, 0x59F111F1B605D019,
    0x923F82A4AF194F9B, 0xAB1C5ED5DA6D8118,
    0xD807AA98A3030242, 0x12835B0145706FBE,
    0x243185BE4EE4B28C, 0x550C7DC3D5FFB4E2,
    0x72BE5D74F27B896F, 0x80DEB1FE3B1696B1,
    0x9BDC06A725C71235, 0xC19BF174CF692694,
    0xE49B69C19EF14AD2, 0xEFBE4786384F25E3,
    0x0FC19DC68B8CD5B5, 0x240CA1CC77AC9C65,
    0x2DE92C6F592B0275, 0x4A7484AA6EA6E483,
    0x5CB0A9DCBD41FBD4, 0x76F988DA831153B5,
    0x983E5152EE66DFAB, 0xA831C66D2DB43210,
    0xB00327C898FB213F, 0xBF597FC7BEEF0EE4,
    0xC6E00BF33DA88FC2, 0xD5A79147930AA725,
    0x06CA6351E003826F, 0x142929670A0E6E70,
    0x27B70A8546D22FFC, 0x2E1B21385C26C926,
    0x4D2C6DFC5AC42AED, 0x53380D139D95B3DF,
    0x650A73548BAF63DE, 0x766A0ABB3C77B2A8,
    0x81C2C92E47EDAEE6, 0x92722C851482353B,
    0xA2BFE8A14CF10364, 0xA81A664BBC423001,
    0xC24B8B70D0F89791, 0xC76C51A30654BE30,
    0xD192E819D6EF5218, 0xD69906245565A910,
    0xF40E35855771202A, 0x106AA07032BBD1B8,
    0x19A4C116B8D2D0C8, 0x1E376C085141AB53,
    0x2748774CDF8EEB99, 0x34B0BCB5E19B48A8,
    0x391C0CB3C5C95A63, 0x4ED8AA4AE3418ACB,
    0x5B9CCA4F7763E373, 0x682E6FF3D6B2B8A3,
    0x748F82EE5DEFB2FC, 0x78A5636F43172F60,
    0x84C87814A1F0AB72, 0x8CC702081A6439EC,
    0x90BEFFFA23631E28, 0xA4506CEBDE82BDE9,
    0xBEF9A3F7B2C67915, 0xC67178F2E372532B,
    0xCA273ECEEA26619C, 0xD186B8C721C0C207,
    0xEADA7DD6CDE0EB1E, 0xF57D4F7FEE6ED178,
    0x06F067AA72176FBA, 0x0A637DC5A2C898A6,
    0x113F9804BEF90DAE, 0x1B710B35131C471B,
    0x28DB77F523047D84, 0x32CAAB7B40C72493,
    0x3C9EBE0A15C9BEBC, 0x431D67C49C100D4C,
    0x4CC5D4BECB3E42B6, 0x597F299CFC657E2A,
    0x5FCB6FAB3AD6FAEC, 0x6C44198C4A475817
};

/*
 * SHA-512 compression function operates on a 1024-bit message block and a
 * 512-bit intermediate hash value (divided in eight 64-bit words).
 * It is essentially a 512-bit block cipher which encrypts the intermediate
 * hash value using the message block as key.
 * (http://www.di-srv.unisa.it/~ads/corso-security/www/CORSO-9900/hash/sha2.pdf)
 */
static void sha512_process(cry_sha512_ctx *ctx, const unsigned char *data)
{
    unsigned int i;
    uint64_t W[80];
    uint64_t a, b, c, d, e, f, g, h;
    uint64_t t1, t2;

    for (i = 0; i < 16; i++)
        CRY_READ64_BE(W[i], data + (i<<3));
    for ( ; i < 80; i++)
        W[i] = S1(W[i - 2]) + W[i - 7] + S0(W[i - 15]) + W[i - 16];

    a = ctx->state[0];
    b = ctx->state[1];
    c = ctx->state[2];
    d = ctx->state[3];
    e = ctx->state[4];
    f = ctx->state[5];
    g = ctx->state[6];
    h = ctx->state[7];

    i = 0;
    do {
        P(a, b, c, d, e, f, g, h, W[i], K[i]);
        i++;
        P(h, a, b, c, d, e, f, g, W[i], K[i]);
        i++;
        P(g, h, a, b, c, d, e, f, W[i], K[i]);
        i++;
        P(f, g, h, a, b, c, d, e, W[i], K[i]);
        i++;
        P(e, f, g, h, a, b, c, d, W[i], K[i]);
        i++;
        P(d, e, f, g, h, a, b, c, W[i], K[i]);
        i++;
        P(c, d, e, f, g, h, a, b, W[i], K[i]);
        i++;
        P(b, c, d, e, f, g, h, a, W[i], K[i]);
        i++;
    } while (i < 80);

    ctx->state[0] += a;
    ctx->state[1] += b;
    ctx->state[2] += c;
    ctx->state[3] += d;
    ctx->state[4] += e;
    ctx->state[5] += f;
    ctx->state[6] += g;
    ctx->state[7] += h;
}


void cry_sha512_init(cry_sha512_ctx *ctx, char is384)
{
    memset(ctx, 0, sizeof(*ctx));

    ctx->is384 = is384;

    /* Load magic initialization constants. */
    if(ctx->is384 == 0) {
        /* SHA-512 */
        ctx->state[0] = 0x6A09E667F3BCC908;
        ctx->state[1] = 0xBB67AE8584CAA73B;
        ctx->state[2] = 0x3C6EF372FE94F82B;
        ctx->state[3] = 0xA54FF53A5F1D36F1;
        ctx->state[4] = 0x510E527FADE682D1;
        ctx->state[5] = 0x9B05688C2B3E6C1F;
        ctx->state[6] = 0x1F83D9ABFB41BD6B;
        ctx->state[7] = 0x5BE0CD19137E2179;
    } else {
        /* SHA-384 */
        ctx->state[0] = 0xCBBB9D5DC1059ED8;
        ctx->state[1] = 0x629A292A367CD507;
        ctx->state[2] = 0x9159015A3070DD17;
        ctx->state[3] = 0x152FECD8F70E5939;
        ctx->state[4] = 0x67332667FFC00B31;
        ctx->state[5] = 0x8EB44A8768581511;
        ctx->state[6] = 0xDB0C2E0D64F98FA7;
        ctx->state[7] = 0x47B5481DBEFA4FA4;
    }
}

void cry_sha512_clear(cry_sha512_ctx *ctx)
{
    cry_memset(ctx, 0, sizeof(*ctx));
}

void cry_sha512_update(cry_sha512_ctx *ctx, const unsigned char *data,
                       size_t len)
{
    size_t fill;
    uint32_t left;

    if (len == 0)
        return;

    left = ctx->len[0] & 0x7F;
    fill = 128 - left;

    ctx->len[0] += (uint32_t)len;
    if (ctx->len[0] < (uint32_t)len)
        ctx->len[1]++;

    if (left && len >= fill) {
        memcpy(ctx->data + left, data, fill);
        sha512_process(ctx, ctx->data);
        data += fill;
        len  -= fill;
        left = 0;
    }

    while (len >= 128) {
        sha512_process(ctx, data);
        data += 128;
        len  -= 128;
    }

    if (len > 0)
        memcpy(ctx->data + left, data, len);
}

void cry_sha512_digest(cry_sha512_ctx *ctx, unsigned char *out)
{
    uint32_t last;
    uint64_t hi, lo;
    size_t padlen;
    unsigned char msglen[16];
    unsigned char pad[128] = { 0 };

    /* Length in bits */
    hi = (ctx->len[0] >> 61) | (ctx->len[1] <<  3);
    lo = (ctx->len[0] <<  3);

    CRY_WRITE64_BE(hi, msglen);
    CRY_WRITE64_BE(lo, msglen + 8);

    last = ctx->len[0] & 0x7F;
    padlen = (last < 112) ? (112 - last) : (128 - last);

    pad[0] = 0x80;
    cry_sha512_update(ctx, pad, padlen);
    cry_sha512_update(ctx, msglen, 16);

    CRY_WRITE64_BE(ctx->state[0], out +  0);
    CRY_WRITE64_BE(ctx->state[1], out +  8);
    CRY_WRITE64_BE(ctx->state[2], out + 16);
    CRY_WRITE64_BE(ctx->state[3], out + 24);
    CRY_WRITE64_BE(ctx->state[4], out + 32);
    CRY_WRITE64_BE(ctx->state[5], out + 40);
    if (ctx->is384 == 0) {
        CRY_WRITE64_BE(ctx->state[6], out + 48);
        CRY_WRITE64_BE(ctx->state[7], out + 56);
    }
}

void cry_sha512(unsigned char *out, const unsigned char *data, size_t len)
{
    cry_sha512_ctx ctx;

    cry_sha512_init(&ctx, 0);
    cry_sha512_update(&ctx, data, len);
    cry_sha512_digest(&ctx, out);
    cry_sha512_clear(&ctx);
}

void cry_sha384(unsigned char *out, const unsigned char *data, size_t len)
{
    cry_sha512_ctx ctx;

    cry_sha512_init(&ctx, 1);
    cry_sha512_update(&ctx, data, len);
    cry_sha512_digest(&ctx, out);
    cry_sha512_clear(&ctx);
}
