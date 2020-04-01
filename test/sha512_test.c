#include "test.h"
#include <cry/sha512.h>

struct sha512_param {
    size_t len;
    unsigned char *data;
    unsigned char hash[CRY_SHA512_DIGEST_SIZE];
};

static void param_init(struct sha512_param *par, int argc, char *argv[])
{
    memset(par, 0, sizeof(*par));
    par->len = strlen(argv[0]) >> 1;
    par->data = malloc(par->len);
    par->len = raw_init(par->data, par->len, argv[0]);
    raw_init(par->hash, CRY_SHA512_DIGEST_SIZE, argv[1]);
}

static void sha384_digest(int argc, char *argv[])
{
    struct sha512_param par;
    unsigned char out[CRY_SHA512_DIGEST_SIZE];

    param_init(&par, argc, argv);
    cry_sha384(out, par.data, par.len);
    ASSERT_EQ_BUF(out, par.hash, 48);
    free(par.data);
}

static void sha512_digest(int argc, char *argv[])
{
    struct sha512_param par;
    unsigned char out[CRY_SHA512_DIGEST_SIZE];

    param_init(&par, argc, argv);
    cry_sha512(out, par.data, par.len);
    ASSERT_EQ_BUF(out, par.hash, CRY_SHA512_DIGEST_SIZE);
    free(par.data);
}

static void dispatch(int argc, char *argv[])
{
    char *test = *argv;

    argv++;
    argc--;
    ASSERT(argc == 2);

    if (strcmp(test, "sha384_digest") == 0)
        sha384_digest(argc, argv);
    else if (strcmp(test, "sha512_digest") == 0)
        sha512_digest(argc, argv);
    else
        TRACE("Test '%s' not defined\n", test);
}

void sha512_test(void)
{
    TRACE("* SHA512 NIST CAVS\n");
    func_test("sha512_test.data", dispatch);
    TRACE("\n");
}
