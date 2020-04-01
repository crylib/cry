#include "test.h"
#include <cry/rsa.h>
#include <cry/prng.h>

static cry_rsa_ctx g_rsa;

static const unsigned char modulus[] = {
    0xC4, 0xF8, 0xE9, 0xE1, 0x5D, 0xCA, 0xDF, 0x2B,
    0x96, 0xC7, 0x63, 0xD9, 0x81, 0x00, 0x6A, 0x64,
    0x4F, 0xFB, 0x44, 0x15, 0x03, 0x0A, 0x16, 0xED,
    0x12, 0x83, 0x88, 0x33, 0x40, 0xF2, 0xAA, 0x0E,
    0x2B, 0xE2, 0xBE, 0x8F, 0xA6, 0x01, 0x50, 0xB9,
    0x04, 0x69, 0x65, 0x83, 0x7C, 0x3E, 0x7D, 0x15,
    0x1B, 0x7D, 0xE2, 0x37, 0xEB, 0xB9, 0x57, 0xC2,
    0x06, 0x63, 0x89, 0x82, 0x50, 0x70, 0x3B, 0x3F
};

static const unsigned char private[] = {
    0x8a, 0x7e, 0x79, 0xf3, 0xfb, 0xfe, 0xa8, 0xeb,
    0xfd, 0x18, 0x35, 0x1c, 0xb9, 0x97, 0x91, 0x36,
    0xf7, 0x05, 0xb4, 0xd9, 0x11, 0x4a, 0x06, 0xd4,
    0xaa, 0x2f, 0xd1, 0x94, 0x38, 0x16, 0x67, 0x7a,
    0x53, 0x74, 0x66, 0x18, 0x46, 0xa3, 0x0c, 0x45,
    0xb3, 0x0a, 0x02, 0x4b, 0x4d, 0x22, 0xb1, 0x5a,
    0xb3, 0x23, 0x62, 0x2b, 0x2d, 0xe4, 0x7b, 0xa2,
    0x91, 0x15, 0xf0, 0x6e, 0xe4, 0x2c, 0x41
};

static const unsigned char public[] = {
    0x01, 0x00, 0x01
};

static const unsigned char plain_text[] = {
    0x09
};

/* Signature for RSA PKCS #1 v15 */
static const unsigned char sign[] = {
    0x26, 0x78, 0xbe, 0x19, 0x97, 0x95, 0x2e, 0x0a,
    0xed, 0xb2, 0xc6, 0xe6, 0xdb, 0xbd, 0x07, 0x7d,
    0x92, 0x1f, 0xca, 0xc0, 0x74, 0x00, 0xe5, 0x6a,
    0x66, 0x8d, 0x95, 0x92, 0xf2, 0xad, 0xc8, 0x3b,
    0x06, 0xda, 0x29, 0x78, 0x1f, 0x48, 0x87, 0xe2,
    0x08, 0x40, 0x03, 0x29, 0x32, 0x78, 0xa9, 0xc6,
    0xee, 0xf9, 0x22, 0x28, 0xce, 0x9e, 0x37, 0x97,
    0x74, 0x46, 0x88, 0x96, 0xf6, 0x24, 0x02, 0xa0
};

/* Signature for Schoolbook RSA */
static const unsigned char sign_sb[] = {
    0xc0, 0x3c, 0x32, 0x14, 0xb1, 0x81, 0x43, 0xbf,
    0xdd, 0xd9, 0x2b, 0x69, 0xe0, 0x1c, 0xd5, 0x19,
    0xa0, 0x06, 0x95, 0xac, 0xdb, 0xc7, 0x16, 0x2d,
    0xd0, 0xac, 0x31, 0x8b, 0x97, 0x78, 0x17, 0x14,
    0x11, 0xf9, 0x57, 0x0b, 0x4a, 0x11, 0x8a, 0x4a,
    0x60, 0xf0, 0xcd, 0xf3, 0xc0, 0x5e, 0xe2, 0xb9,
    0x5c, 0x47, 0x95, 0xa4, 0xd8, 0x73, 0xc2, 0x29,
    0x9d, 0x58, 0x77, 0x47, 0xf3, 0x05, 0x78, 0x79
};

#define PLAIN_LEN sizeof(plain_text)
#define SIGN_LEN sizeof(sign)
#define CIPHER_LEN SIGN_LEN

static void setup(void)
{
    /* Seed the PRNG to make results predictable */
    cry_prng_aes_init(NULL, 0);
    cry_rsa_init(&g_rsa, CRY_RSA_PADDING_PKCS_V15);
}

static void teardown(void)
{
    cry_rsa_clear(&g_rsa);
}

static void encrypt_decrypt(void)
{
    size_t outlen;
    unsigned char *cipher_buf;
    unsigned char *plain_buf;

    cry_mpi_load_bin(&g_rsa.n, modulus, sizeof(modulus));
    cry_mpi_load_bin(&g_rsa.e, public, sizeof(public));
    cry_mpi_load_bin(&g_rsa.d, private, sizeof(private));

    ASSERT_OK(cry_rsa_encrypt(&g_rsa, &cipher_buf, &outlen,
                              plain_text, PLAIN_LEN));
    ASSERT_EQ(outlen, CIPHER_LEN);

    ASSERT_OK(cry_rsa_decrypt(&g_rsa, &plain_buf, &outlen,
                              cipher_buf, outlen));
    ASSERT_EQ(outlen, PLAIN_LEN);
    ASSERT_EQ_BUF(plain_buf, plain_text, outlen);
    free(plain_buf);
    free(cipher_buf);
}

static void encrypt_decrypt_sb(void)
{
    size_t outlen, plain_len;
    unsigned char *cipher_buf;
    unsigned char *plain_buf;
    unsigned char *plain_text;

    cry_mpi_load_bin(&g_rsa.n, modulus, sizeof(modulus));
    cry_mpi_load_bin(&g_rsa.e, public, sizeof(public));
    cry_mpi_load_bin(&g_rsa.d, private, sizeof(private));

    g_rsa.padding = CRY_RSA_PADDING_NONE;

    plain_len = 10;
    plain_text = malloc(plain_len);
    memset(plain_text, 'A', plain_len);
    ASSERT_OK(cry_rsa_encrypt(&g_rsa, &cipher_buf, &outlen,
                              plain_text, plain_len));
    ASSERT_EQ(outlen, CIPHER_LEN);

    ASSERT_OK(cry_rsa_decrypt(&g_rsa, &plain_buf, &outlen,
                              cipher_buf, outlen));
    ASSERT_EQ(outlen, cry_mpi_count_bytes(&g_rsa.n));
    ASSERT_EQ_BUF(plain_buf, plain_text, plain_len);
    free(plain_buf);
    free(cipher_buf);
    free(plain_text);
}

static void sign_verify(void)
{
    size_t siglen;
    unsigned char *sig_buf;

    cry_mpi_load_bin(&g_rsa.n, modulus, sizeof(modulus));
    cry_mpi_load_bin(&g_rsa.e, public, sizeof(public));
    cry_mpi_load_bin(&g_rsa.d, private, sizeof(private));

    ASSERT_OK(cry_rsa_sign(&g_rsa, &sig_buf, &siglen,
                           plain_text, PLAIN_LEN));
    ASSERT_EQ(siglen, CIPHER_LEN);
    ASSERT_EQ_BUF(sig_buf, sign, siglen);

    ASSERT_OK(cry_rsa_verify(&g_rsa, sig_buf, siglen,
                             plain_text, PLAIN_LEN));
    free(sig_buf);
}

static void sign_verify_sb(void)
{
    size_t siglen;
    unsigned char *sig_buf;

    cry_mpi_load_bin(&g_rsa.n, modulus, sizeof(modulus));
    cry_mpi_load_bin(&g_rsa.e, public, sizeof(public));
    cry_mpi_load_bin(&g_rsa.d, private, sizeof(private));

    g_rsa.padding = CRY_RSA_PADDING_NONE;

    ASSERT_OK(cry_rsa_sign(&g_rsa, &sig_buf, &siglen,
                           plain_text, PLAIN_LEN));
    ASSERT_EQ(siglen, CIPHER_LEN);
    ASSERT_EQ_BUF(sig_buf, sign_sb, siglen);

    ASSERT_OK(cry_rsa_verify(&g_rsa, sig_buf, siglen,
                             plain_text, PLAIN_LEN));
    free(sig_buf);
}

struct rsa_param {
    size_t mlen;
    size_t elen;
    size_t clrlen;
    size_t ciplen;
    unsigned char *mraw;
    unsigned char *eraw;
    unsigned char *clrraw;
    unsigned char *cipraw;
};

/*
 * Params:
 * p0 : M
 * p1 : E
 * p2 : cleartext
 * p3 : ciphertext
 */
static void rsa_param_init(struct rsa_param *par, int argc, char *argv[])
{
    memset(par, 0, sizeof(*par));

    ASSERT(argc == 4);

    par->mlen = strlen(argv[0]) >> 1;
    par->elen = strlen(argv[1]) >> 1;
    par->clrlen = strlen(argv[2]) >> 1;
    par->ciplen = strlen(argv[3]) >> 1;
    par->mraw = malloc(par->mlen + par->elen + par->clrlen + par->ciplen);
    par->eraw = par->mraw + par->mlen;
    par->clrraw = par->eraw + par->elen;
    par->cipraw = par->clrraw + par->clrlen;

    raw_init(par->mraw, par->mlen, argv[0]);
    raw_init(par->eraw, par->elen, argv[1]);
    raw_init(par->clrraw, par->clrlen, argv[2]);
    raw_init(par->cipraw, par->ciplen, argv[3]);
}

static void rsa_pkcs1_encrypt(int argc, char *argv[])
{
    struct rsa_param par;
    size_t outlen;
    unsigned char *cipher_buf;

    rsa_param_init(&par, argc, argv);
    cry_mpi_load_bin(&g_rsa.n, par.mraw, par.mlen);
    cry_mpi_load_bin(&g_rsa.e, par.eraw, par.elen);

    ASSERT_OK(cry_rsa_encrypt(&g_rsa, &cipher_buf, &outlen,
                              par.clrraw, par.clrlen));

    free(par.mraw);
    free(cipher_buf);
}

static void dispatch(int argc, char *argv[])
{
    char *test = *argv;

    setup();
    argv++;
    argc--;
    if (strcmp(test, "rsa_pkcs1_encrypt") == 0)
        rsa_pkcs1_encrypt(argc, argv);
    else
        TRACE("Test '%s' not defined\n", test);
    teardown();
}

/* Just for... coverage :-) */
#define KEYGEN_BITS 256

static void keygen_known_exp(void)
{
    ASSERT_OK(cry_rsa_keygen(&g_rsa, KEYGEN_BITS, 65537));
}

static void keygen_rand_exp(void)
{
    ASSERT_OK(cry_rsa_keygen(&g_rsa, KEYGEN_BITS, 0));
}

#define MYRUN(name, test) \
        run(name, test, setup, teardown)

void rsa_test(void)
{
    TRACE("* RSA Keygen\n");
    MYRUN("Keygen 512 with 65537 exponent", keygen_known_exp);
    MYRUN("Keygen 512 with random exponent", keygen_rand_exp);

    TRACE("* RSA Encryption\n");
    MYRUN("Encrypt-Decrypt Schoolbook", encrypt_decrypt_sb);
    MYRUN("Encrypt-Decrypt PKCS v15", encrypt_decrypt);

    TRACE("* RSA Digital Signature\n");
    MYRUN("Sign-Verify Schoolbook", sign_verify_sb);
    MYRUN("Sign-Verify PKCS v15", sign_verify);
    func_test("rsa_test.data", dispatch);

    TRACE("\n");
}
