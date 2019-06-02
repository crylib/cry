#include <cry/sha256.h>
#include <string.h>
#include "../misc.h"

#define ROTLEFT(a, b)  (((a) << (b)) | ((a) >> (32 - (b))))
#define ROTRIGHT(a, b) (((a) >> (b)) | ((a) << (32 - (b))))

#define CH(x, y, z)  (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x, y, z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x)  (ROTRIGHT(x,  2) ^ ROTRIGHT(x, 13) ^ ROTRIGHT(x, 22))
#define EP1(x)  (ROTRIGHT(x,  6) ^ ROTRIGHT(x, 11) ^ ROTRIGHT(x, 25))
#define SIG0(x) (ROTRIGHT(x,  7) ^ ROTRIGHT(x, 18) ^ ((x) >> 3))
#define SIG1(x) (ROTRIGHT(x, 17) ^ ROTRIGHT(x, 19) ^ ((x) >> 10))

static const uint32_t k[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

static void cry_sha256_transform(cry_sha256_ctx *ctx,
                                 const unsigned char *data)
{
    uint32_t a, b, c, d, e, f, g, h, t1, t2, m[64];
    unsigned int i, j;

    for (i = 0, j = 0; i < 16; ++i, j += 4)
        CRY_READ32_BE(m[i], &data[j]);

    for (; i < 64; ++i)
        m[i] = SIG1(m[i - 2]) + m[i - 7] + SIG0(m[i - 15]) + m[i - 16];

    a = ctx->state[0];
    b = ctx->state[1];
    c = ctx->state[2];
    d = ctx->state[3];
    e = ctx->state[4];
    f = ctx->state[5];
    g = ctx->state[6];
    h = ctx->state[7];

    for (i = 0; i < 64; ++i) {
        t1 = h + EP1(e) + CH(e, f, g) + k[i] + m[i];
        t2 = EP0(a) + MAJ(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }

    ctx->state[0] += a;
    ctx->state[1] += b;
    ctx->state[2] += c;
    ctx->state[3] += d;
    ctx->state[4] += e;
    ctx->state[5] += f;
    ctx->state[6] += g;
    ctx->state[7] += h;
}

void cry_sha256_init(cry_sha256_ctx *ctx)
{
    memset(ctx, 0, sizeof(*ctx));
    /* Load magic initialization constants. */
    ctx->state[0] = 0x6a09e667;
    ctx->state[1] = 0xbb67ae85;
    ctx->state[2] = 0x3c6ef372;
    ctx->state[3] = 0xa54ff53a;
    ctx->state[4] = 0x510e527f;
    ctx->state[5] = 0x9b05688c;
    ctx->state[6] = 0x1f83d9ab;
    ctx->state[7] = 0x5be0cd19;
}

void cry_sha256_clear(cry_sha256_ctx *ctx)
{
    cry_memset(ctx, 0, sizeof(*ctx));
}

void cry_sha256_update(cry_sha256_ctx *ctx, const unsigned char *data,
                       size_t size)
{
    size_t i;

    for (i = 0; i < size; ++i) {
        ctx->data[ctx->datalen] = data[i];
        ctx->datalen++;
        if (ctx->datalen == 64) {
            cry_sha256_transform(ctx, ctx->data);
            ctx->bitlen += 512;
            ctx->datalen = 0;
        }
    }
}

void cry_sha256_digest(cry_sha256_ctx *ctx, unsigned char *digest)
{
    size_t i = ctx->datalen;

    /* Pad whatever data is left in the buffer */
    if (ctx->datalen < 56) {
        ctx->data[i++] = 0x80;
        while (i < 56)
            ctx->data[i++] = 0x00;
    } else {
        ctx->data[i++] = 0x80;
        while (i < 64)
            ctx->data[i++] = 0x00;
        cry_sha256_transform(ctx, ctx->data);
        memset(ctx->data, 0, 56);
    }

    /*
     * Append to the padding the total message's length in bits
     * and transform.
     */
    ctx->bitlen += ctx->datalen * 8;
    CRY_WRITE64_BE(ctx->bitlen, &ctx->data[56]);
    cry_sha256_transform(ctx, ctx->data);

    /*
     * Since this implementation uses little endian byte ordering and
     * SHA uses big endian, reverse all the bytes when copying the final
     * state to the output hash.
     */
    for (i = 0; i < 4; ++i) {
        digest[i]      = (ctx->state[0] >> (24 - i * 8)) & 0xff;
        digest[i + 4]  = (ctx->state[1] >> (24 - i * 8)) & 0xff;
        digest[i + 8]  = (ctx->state[2] >> (24 - i * 8)) & 0xff;
        digest[i + 12] = (ctx->state[3] >> (24 - i * 8)) & 0xff;
        digest[i + 16] = (ctx->state[4] >> (24 - i * 8)) & 0xff;
        digest[i + 20] = (ctx->state[5] >> (24 - i * 8)) & 0xff;
        digest[i + 24] = (ctx->state[6] >> (24 - i * 8)) & 0xff;
        digest[i + 28] = (ctx->state[7] >> (24 - i * 8)) & 0xff;
    }
}

void cry_sha256(unsigned char *out, const unsigned char *data, size_t len)
{
    cry_sha256_ctx ctx;

    cry_sha256_init(&ctx);
    cry_sha256_update(&ctx, data, len);
    cry_sha256_digest(&ctx, out);
    cry_sha256_clear(&ctx);
}
