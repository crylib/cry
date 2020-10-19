/**
 * @file    version.h
 *
 * @brief   Software version information.
 *
 * CRY follows the Semantic Versioning 2.0.0 (https://semver.org).
 *
 * Given a version number MAJOR.MINOR.PATCH, increment the:
 * - MAJOR version when you make incompatible API changes,
 * - MINOR version when you add functionality in a backwards compatible manner,
 * - PATCH version when you make backwards compatible bug fixes.
 * - Additional labels for pre-release and build metadata are available as
 *   extensions to the MAJOR.MINOR.PATCH format.
 */

#ifndef CRY_VERSION_H_
#define CRY_VERSION_H_

/** Library pre-release string */
#define CRY_PRE_RELEASE "alpha.1"

/** Library version major number */
#define CRY_MAJOR   0
/** Library version minor number */
#define CRY_MINOR   1
/** Library version patch number */
#define CRY_PATCH   2

/** Stable version */
#ifdef CRY_PRE_RELEASE
#define CRY_STABLE  0
#else
#define CRY_STABLE  1
#endif

/** Library version number. */
#define CRY_VERSION \
    ((long)CRY_MAJOR << 28 | \
     (long)CRY_MINOR << 20 | \
     (long)CRY_PATCH <<  4 | \
     (long)CRY_STABLE)

#define _CRY_STR_HELPER(x) #x
#define _CRY_STR(x) _CRY_STR_HELPER(x)
#define _CRY_VERSION_STR \
    _CRY_STR(CRY_MAJOR) "." \
    _CRY_STR(CRY_MINOR) "." \
    _CRY_STR(CRY_PATCH)

#ifdef CRY_PRE_RELEASE
#define CRY_VERSION_STR _CRY_VERSION_STR "-" CRY_PRE_RELEASE
#else
#define CRY_VERSION_STR _CRY_VERSION_STR
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Get library version number.
 *
 * @return  Library version number (CRY_VERSION) at build time.
 */
long cry_version(void);

/**
 * Get library version string.
 *
 * @return  Library version string (CRY_VERSION_STR) at build time.
 */
const char *cry_version_str(void);

#ifdef __cplusplus
}
#endif

#endif /* CRY_VERSION_H_ */
