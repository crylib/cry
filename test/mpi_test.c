#include "test.h"
#include "malloc_mock.h"
#include <cry/mpi.h>
#include <stdlib.h>

#define ERROR_FLAG      'Z'
#define MPI_BUF_LEN     6

static int g_exp_res;
static cry_mpi *g_mpi_buf;

#define g_mpi0 (&g_mpi_buf[0])
#define g_mpi1 (&g_mpi_buf[1])
#define g_mpi2 (&g_mpi_buf[2])
#define g_mpi3 (&g_mpi_buf[3])
#define g_mpi_r1 (&g_mpi_buf[4])
#define g_mpi_r2 (&g_mpi_buf[5])


static void mpi_setup(void)
{
    g_exp_res = 0;
    g_mpi_buf = malloc(sizeof(cry_mpi)*MPI_BUF_LEN);
    cry_mpi_init_list(g_mpi0, g_mpi1, g_mpi2, g_mpi3, g_mpi_r1, g_mpi_r2, NULL);
}

static void mpi_teardown(void)
{
    cry_mpi_clear_list(g_mpi0, g_mpi1, g_mpi2, g_mpi3, g_mpi_r1, g_mpi_r2,
                       NULL);
    free(g_mpi_buf);
}

static int mpi_load_args(int argc, char *argv[])
{
    size_t i;
    int res = 0;

    if (*argv[argc-1] == ERROR_FLAG) {
        argc--;
        g_exp_res = atoi(argv[argc]+1);
    }

    for (i = 0; i < argc; i++) {
        if ((res = cry_mpi_load_str(&g_mpi_buf[i], 16, argv[i])) < 0) {
            while (i-- > 0)
                cry_mpi_clear(&g_mpi_buf[i]);
            break;
        }
    }
    return res;
}

static void mpi_init_size(int argc, char *argv[])
{
    cry_mpi tmp;

    ASSERT_EQ(argc, 1);
    ASSERT_OK(cry_mpi_init_size(&tmp, atol(argv[0])));

    ASSERT_EQ(tmp.alloc, atol(argv[0]));

    cry_mpi_clear(&tmp);
}

static void mpi_count_bits(int argc, char *argv[])
{
    ASSERT_EQ(argc, 2);
    ASSERT_OK(cry_mpi_load_str(g_mpi0, 16, argv[0]));

    ASSERT_EQ(cry_mpi_count_bits(g_mpi0), atol(argv[1]));
}

static void mpi_count_bytes(int argc, char *argv[])
{
    ASSERT_EQ(argc, 2);
    ASSERT_OK(cry_mpi_load_str(g_mpi0, 16, argv[0]));

    ASSERT_EQ(cry_mpi_count_bytes(g_mpi0), atol(argv[1]));
}

static void mpi_set_int(int argc, char *argv[])
{
    long l1, l2;

    ASSERT_EQ(argc, 1);
    l1 = atol(argv[0]);

    ASSERT_OK(cry_mpi_set_int(g_mpi0, l1));
    ASSERT_OK(cry_mpi_get_int(g_mpi0, &l2));

    ASSERT_EQ(l1, l2);
}

static void mpi_load_store_str(int argc, char *argv[])
{
    int res;

    ASSERT_EQ(argc, 3);
    res = atoi(argv[2]);

    ASSERT(cry_mpi_load_str(g_mpi0, 16, argv[0]) == res);
    if (res == 0) {
        ASSERT(cry_mpi_store_str(g_mpi0, 16, (char *)g_buf) == 0);
        ASSERT(strcmp((char *)g_buf, argv[1]) == 0);
    }
}

static void mpi_cmp(int argc, char *argv[])
{
    int res;

    ASSERT_EQ(argc, 3);
    res = atoi(argv[2]);

    ASSERT(cry_mpi_load_str(g_mpi0, 16, argv[0]) == 0);
    ASSERT(cry_mpi_load_str(g_mpi1, 16, argv[1]) == 0);

    ASSERT(cry_mpi_cmp(g_mpi0, g_mpi1) == res);
}


static void mpi_shift(int argc, char *argv[], int left)
{
    int res, bits;

    ASSERT_EQ(argc, 3);
    bits = atoi(argv[1]);
    ASSERT(cry_mpi_load_str(g_mpi0, 16, argv[0]) == 0);
    ASSERT(cry_mpi_load_str(g_mpi1, 16, argv[2]) == 0);

    if (left == 1)
        res = cry_mpi_shl(g_mpi_r1, g_mpi0, bits);
    else
        res = cry_mpi_shr(g_mpi_r1, g_mpi0, bits);

    ASSERT(res == 0);
    ASSERT(cry_mpi_cmp(g_mpi_r1, g_mpi1) == 0);
}


typedef int (*unary_op_f)(cry_mpi *r, const cry_mpi *a);

static void mpi_unary_op(int argc, char *argv[], unary_op_f op)
{
    int res;

    ASSERT_EQ(argc, 2);
    ASSERT_OK(mpi_load_args(argc, argv));

    res = op(g_mpi_r1, g_mpi0);

    ASSERT_EQ(res, g_exp_res);
    if (res == 0)
        ASSERT(cry_mpi_cmp(g_mpi_r1, g_mpi1) == 0);
}


#define RUNM(test) do { \
    if (g_malloc_mock_state == MALLOC_MOCK_READY) \
        g_malloc_mock_state = MALLOC_MOCK_ACTIVE; \
    test; \
    if (g_malloc_mock_state == MALLOC_MOCK_FAILED) { \
        ASSERT(res == -1); \
        return; \
    } \
} while (0)


typedef int (*binary_op_f)(cry_mpi *r, const cry_mpi *a, const cry_mpi *b);

static void mpi_binary_op(int argc, char *argv[], binary_op_f op)
{
    int res;

    ASSERT_EQ(argc, 3);
    ASSERT_OK(mpi_load_args(argc, argv));

    RUNM(res = op(g_mpi_r1, g_mpi0, g_mpi1));

    ASSERT_EQ(res, g_exp_res);
    if (res == 0)
        ASSERT(cry_mpi_cmp(g_mpi_r1, g_mpi2) == 0);
}


typedef int (*binary_mod_op_f)(cry_mpi *r, const cry_mpi *a,
                               const cry_mpi *b, const cry_mpi *m);

static void mpi_binary_mod_op(int argc, char *argv[], binary_mod_op_f op)
{
    int res;

    ASSERT_EQ(argc, 4);
    ASSERT_OK(mpi_load_args(argc, argv));

    RUNM(res = op(g_mpi_r1, g_mpi0, g_mpi1, g_mpi2));

    ASSERT_EQ(res, g_exp_res);
    if (res == 0)
        ASSERT(cry_mpi_cmp(g_mpi_r1, g_mpi3) == 0);
}


static void mpi_div(int argc, char *argv[])
{
    int res;

    ASSERT_EQ(argc, 4);
    ASSERT(mpi_load_args(argc, argv) == 0);

    RUNM(res = cry_mpi_div(g_mpi_r1, g_mpi_r2, g_mpi0, g_mpi1));

    ASSERT(res == g_exp_res);
    if (res == 0) {
        ASSERT(cry_mpi_cmp(g_mpi_r1, g_mpi2) == 0);
        ASSERT(cry_mpi_cmp(g_mpi_r2, g_mpi3) == 0);
    }
}


static void mpi_dispatch(int argc, char *argv[])
{
    char *test = *argv;

    mpi_setup();

    argv++;
    argc--;
    if (strcmp(test, "mpi_init_size") == 0)
        mpi_init_size(argc, argv);
    else if (strcmp(test, "mpi_count_bits") == 0)
        mpi_count_bits(argc, argv);
    else if (strcmp(test, "mpi_count_bytes") == 0)
        mpi_count_bytes(argc, argv);
    else if (strcmp(test, "mpi_set_int") == 0)
        mpi_set_int(argc, argv);
    else if (strcmp(test, "mpi_load_store_str") == 0)
        mpi_load_store_str(argc, argv);
    else if (strcmp(test, "mpi_abs") == 0)
        mpi_unary_op(argc, argv, cry_mpi_abs);
    else if (strcmp(test, "mpi_cmp") == 0)
        mpi_cmp(argc, argv);
    else if (strcmp(test, "mpi_shl") == 0)
        mpi_shift(argc, argv, 1);
    else if (strcmp(test, "mpi_shr") == 0)
        mpi_shift(argc, argv, 0);
    else if (strcmp(test, "mpi_add") == 0)
        mpi_binary_op(argc, argv, cry_mpi_add);
    else if (strcmp(test, "mpi_sub") == 0)
        mpi_binary_op(argc, argv, cry_mpi_sub);
    else if (strcmp(test, "mpi_mul") == 0)
        mpi_binary_op(argc, argv, cry_mpi_mul);
    else if (strcmp(test, "mpi_mul_baseline") == 0)
        mpi_binary_op(argc, argv, cry_mpi_mul_baseline);
    else if (strcmp(test, "mpi_mul_comba") == 0)
        mpi_binary_op(argc, argv, cry_mpi_mul_comba);
    else if (strcmp(test, "mpi_mul_karatsuba") == 0)
        mpi_binary_op(argc, argv, cry_mpi_mul_karatsuba);
    else if (strcmp(test, "mpi_mul_toom3") == 0)
        mpi_binary_op(argc, argv, cry_mpi_mul_toom3);
    else if (strcmp(test, "mpi_div") == 0)
        mpi_div(argc, argv);
    else if (strcmp(test, "mpi_exp") == 0)
        mpi_binary_op(argc, argv, cry_mpi_exp);
    else if (strcmp(test, "mpi_mod_exp") == 0)
        mpi_binary_mod_op(argc, argv, cry_mpi_mod_exp);
    else if (strcmp(test, "mpi_sqr") == 0)
        mpi_unary_op(argc, argv, cry_mpi_sqr);
    else if (strcmp(test, "mpi_sqr_baseline") == 0)
        mpi_unary_op(argc, argv, cry_mpi_sqr_baseline);
    else if (strcmp(test, "mpi_sqr_comba") == 0)
        mpi_unary_op(argc, argv, cry_mpi_sqr_comba);
    else if (strcmp(test, "mpi_sqr_karatsuba") == 0)
        mpi_unary_op(argc, argv, cry_mpi_sqr_karatsuba);
    else
        TRACE("Test '%s' not defined\n", test);

    mpi_teardown();
}


static struct malloc_fail_args g_malloc_fail_tests[] = {
    { 5,
      { "Multiply Baseline", "mpi_mul_baseline", "12345678", "12345678",
        "14b66dc1df4d840" }},
    { 5,
      { "Multiply Comba", "mpi_mul_comba", "12345678", "12345678",
        "14b66dc1df4d840" }},
    { 5,
      { "Multiply Karatsuba", "mpi_mul_karatsuba", "FFFFFFFFFFFFFFFF",
        "FFFFFFFFFFFFFFFF", "fffffffffffffffe0000000000000001" }},
    { 5,
      { "Multiply Toom-3", "mpi_mul_toom3", "FFFFFFFFFFFFFFFF",
        "FFFFFFFFFFFFFFFF", "fffffffffffffffe0000000000000001" }},
};


void mpi_test(void)
{
    TRACE("* MPI\n");
    func_test("mpi_test.data", mpi_dispatch);
    malloc_fail_tests(g_malloc_fail_tests,
                      ARLEN(g_malloc_fail_tests), mpi_dispatch);
    TRACE("\n");
}
