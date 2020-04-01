#include <cry/dh.h>
#include "misc.h"


int cry_dh_init(cry_dh_ctx *ctx)
{
    return cry_mpi_init_list(&ctx->p, &ctx->g, &ctx->e, &ctx->Y, &ctx->Z,
                             (cry_mpi *)NULL);
}

void cry_dh_clear(cry_dh_ctx *ctx)
{
    cry_mpi_clear_list(&ctx->p, &ctx->g, &ctx->e, &ctx->Y, &ctx->Z,
                       (cry_mpi *)NULL);
    cry_memset(ctx, 0, sizeof(*ctx));
}

int cry_dh_agree(cry_dh_ctx *ctx)
{
    return cry_mpi_mod_exp(&ctx->Y, &ctx->g, &ctx->e, &ctx->p);
}

int cry_dh_finalize(cry_dh_ctx *ctx)
{
    return cry_mpi_mod_exp(&ctx->Z, &ctx->Y, &ctx->e, &ctx->p);
}

int cry_dh_get_tok(cry_dh_ctx *ctx, unsigned char *out, size_t out_len)
{
    return cry_mpi_store_bin(&ctx->Y, out, out_len, 1);
}

int cry_dh_set_tok(cry_dh_ctx *ctx, unsigned char *in, size_t in_len)
{
    return cry_mpi_load_bin(&ctx->Y, in, in_len);
}

int cry_dh_get_sec(cry_dh_ctx *ctx, unsigned char *out, size_t out_len)
{
    return cry_mpi_store_bin(&ctx->Z, out, out_len, 1);
}
