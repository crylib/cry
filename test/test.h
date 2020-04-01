#ifndef CRY_TEST_H_
#define CRY_TEST_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cry/config.h>

#define BIGBUF_SIZ  8192

/* Constant seed to make "rand" results predictable */
#define RAND_SEED_RAW  ((unsigned char *)"RANDSEED")
#define RAND_SEED_SIZ  (sizeof(RAND_SEED_RAW) - 1)

#define ARLEN(ar) (sizeof(ar)/sizeof(ar[0]))

/* Skip tests that are extremely slow when using 1 byte digits */
#ifdef CRY_MPI_DEBUG_CONF
#define SKIP_SLOW
#endif

extern int g_fails;
extern unsigned char g_buf[BIGBUF_SIZ];
extern int g_verbose;

typedef void (*test_func_t)(void);

struct test_case {
    const char *name;
    test_func_t func;
};


typedef void (*dispatch_func_t)(int argc, char *argv[]);

void func_test(const char *datafile, dispatch_func_t dispatch);


struct malloc_fail_args {
    int argc;
    char *argv[10];
};

void malloc_fail_tests(struct malloc_fail_args *args, size_t num,
                       dispatch_func_t dispatch);

void asc_to_raw(const char *asc, size_t size, unsigned char *raw);

int raw_init(unsigned char *raw, size_t rawlen, const char *asc);


void run(const char *name, void (*test)(void),
         void (*setup)(void), void (*teardown)(void));

#define RUN(test) run(#test, test, NULL, NULL)


#define TRACE(...) do { \
    fprintf(stdout, __VA_ARGS__); \
    fflush(stdout); \
} while(0)

#define TRACE2(...) do { \
    if (g_verbose) \
        TRACE(__VA_ARGS__); \
} while(0)

#define PRINT_HEX(msg, buf, siz) do { \
    size_t __siz = siz; \
    unsigned char *__p = buf; \
    if (msg) \
        TRACE("%s: ", msg); \
    while (__siz-- > 0) \
        TRACE("%02x", *__p++); \
    TRACE("\n"); \
} while(0)

#define PRINT_ASC(msg, buf, siz) \
    TRACE("%s: %.*s\n", msg, (int)(siz), buf)

#define PRINT_MPI(msg, mpi, rad) do { \
    TRACE("%s:\t", msg); \
    cry_mpi_print(mpi, rad); \
} while(0)


#define ASSERT(test) do { \
    if ((test) == 0) { \
        g_fails++; \
        TRACE("!!! ASSERTION FAIL (%s:%d): %s\n", \
                __FILE__, __LINE__, #test); \
        return; \
    } \
    } while (0)

#define ASSERT_EQ(i1, i2) \
    ASSERT((i1) == (i2))

#define ASSERT_NE(i1, i2) \
    ASSERT((i1) != (i2))

#define ASSERT_EQ_BUF(b1, b2, len) \
    ASSERT(memcmp(b1, b2, len) == 0)

#define ASSERT_NE_BUF(b1, b2, len) \
    ASSERT(memcmp(b1, b2, len) != 0)

#define ASSERT_OK(e) \
    ASSERT_EQ((e), 0)


#endif /* _TEST_H_ */
