#include "test.h"
#include <cry/version.h>

void version_test_num(void)
{
    ASSERT_EQ(cry_version(), CRY_VERSION);
}

void version_test_str(void)
{
    ASSERT_EQ(strcmp(cry_version_str(), CRY_VERSION_STR), 0);
}

void version_test(void)
{
    TRACE("* Version\n");
    RUN(version_test_num);
    RUN(version_test_str);
    TRACE("\n");
}
