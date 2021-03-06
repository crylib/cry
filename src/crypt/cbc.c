#include <cry/cbc.h>
#include <cry/ciph.h>
#include "memxor.h"
#include "misc.h"
#include <string.h>
#include <stdlib.h>


void cry_cbc_init(struct cry_cbc_ctx *ctx, void *ciph_ctx,
                  const cry_ciph_itf *ciph_itf)
{
    memset(ctx, 0, sizeof(*ctx));
    ctx->ciph_ctx = ciph_ctx;
    ctx->ciph_itf = ciph_itf;
    if (ciph_itf->init != NULL)
        ciph_itf->init(ciph_ctx);
}

void cry_cbc_clear(struct cry_cbc_ctx *ctx)
{
    if (ctx->ciph_itf->clear != NULL)
        ctx->ciph_itf->clear(ctx->ciph_ctx);
    memset(ctx, 0, sizeof(*ctx));
}

void cry_cbc_key_set(struct cry_cbc_ctx *ctx, const unsigned char *key,
                     size_t size)
{
    void *ciph = ctx->ciph_ctx;
    cry_ciph_key_set_f key_set = ctx->ciph_itf->key_set;

    key_set(ciph, key, size);
}

void cry_cbc_iv_set(struct cry_cbc_ctx *ctx, const unsigned char *iv,
                    size_t size)
{
    ctx->ctrlen = CRY_MIN(CRY_CBC_BLOCK_MAX, size);
    memcpy(ctx->ctr, iv, ctx->ctrlen);
    memset(ctx->ctr + ctx->ctrlen, 0, CRY_CBC_BLOCK_MAX - ctx->ctrlen);
}

void cry_cbc_encrypt(struct cry_cbc_ctx *ctx, unsigned char *dst,
                     const unsigned char *src, size_t size)
{
    void *ciph = ctx->ciph_ctx;
    cry_ciph_encrypt_f encrypt = ctx->ciph_itf->encrypt;

    size &= ~(ctx->ctrlen - 1);

    while (size != 0) {
        cry_memxor(ctx->ctr, src, ctx->ctrlen);
        encrypt(ciph, dst, ctx->ctr, ctx->ctrlen);
        memcpy(ctx->ctr, dst, ctx->ctrlen);

        src += ctx->ctrlen;
        dst += ctx->ctrlen;
        size -= ctx->ctrlen;
    }
}

#define CBC_BUFFER_LIMIT  512

void cry_cbc_decrypt(struct cry_cbc_ctx *ctx, unsigned char *dst,
                     const unsigned char *src, size_t size)
{
    void *ciph = ctx->ciph_ctx;
    cry_ciph_encrypt_f decrypt = ctx->ciph_itf->decrypt;

    size &= ~(ctx->ctrlen - 1);

    if (src != dst) {
        decrypt(ciph, dst, src, size);
        cry_memxor(dst, ctx->ctr, ctx->ctrlen);
        cry_memxor(dst + ctx->ctrlen, src, size - ctx->ctrlen);
        memcpy(ctx->ctr, src + size - ctx->ctrlen, ctx->ctrlen);
    } else {
        /*
         * For in-place CBC, we decrypt into a temporary buffer of size
         * at most CBC_BUFFER_LIMIT, and process that amount of data at
         * a time.
         */
        unsigned char init_iv[CRY_CBC_BLOCK_MAX];
        unsigned char fallback[CRY_CBC_BLOCK_MAX]; /* if malloc fails */
        unsigned char *buf;
        size_t buf_siz = CRY_MIN(size, CBC_BUFFER_LIMIT);

        buf = malloc(buf_siz);
        if (buf == NULL) {
            buf = fallback;
            buf_siz = CRY_CBC_BLOCK_MAX;
        }

        while (size > buf_siz) {
            decrypt(ciph, buf, src, buf_siz);
            memcpy(init_iv, ctx->ctr, ctx->ctrlen);
            memcpy(ctx->ctr, src + buf_siz - ctx->ctrlen, ctx->ctrlen);
            cry_memxor2(dst + ctx->ctrlen, buf + ctx->ctrlen,
                        src, buf_siz - ctx->ctrlen);
            cry_memxor2(dst, buf, init_iv, ctx->ctrlen);
            size -= buf_siz;
            src  += buf_siz;
            dst  += buf_siz;
        }

        decrypt(ciph, buf, src, size);
        memcpy(init_iv, ctx->ctr, ctx->ctrlen);
        /* Copies last block */
        memcpy(ctx->ctr, src + size - ctx->ctrlen, ctx->ctrlen);
        /* Writes all but first block, reads all but last block. */
        cry_memxor2(dst + ctx->ctrlen, buf + ctx->ctrlen,
                    src, size - ctx->ctrlen);
        /* Writes first block. */
        cry_memxor2(dst, buf, init_iv, ctx->ctrlen);

        if (buf != fallback)
            free(buf);
    }
}
