#include "test.h"
#include <cry/base64.h>


static void encode(int argc, char *argv[])
{
    int len;

    ASSERT(argc == 3);
    len = atoi(argv[2]);

    ASSERT(cry_base64_encode((char *)g_buf, argv[0], strlen(argv[0])) == len);
    if (len >= 0)
        ASSERT_EQ_BUF(g_buf, argv[1], len);
}

static void decode(int argc, char *argv[])
{
    int len;

    ASSERT(argc == 3);
    len = atoi(argv[2]);

    ASSERT(cry_base64_decode((char *)g_buf, argv[0], strlen(argv[0])) == len);
    if (len >= 0)
        ASSERT_EQ_BUF(g_buf, argv[1], len);
}


static void base64_dispatch(int argc, char *argv[])
{
    char *test = *argv;

    argv++;
    argc--;
    if (strcmp(test, "encode") == 0)
        encode(argc, argv);
    else if (strcmp(test, "decode") == 0)
        decode(argc, argv);
    else
        TRACE("Test '%s' not defined\n", test);
}


void base64_test(void)
{
    TRACE("* Base-64\n");
    func_test("base64_test.data", base64_dispatch);
    TRACE("\n");
}
