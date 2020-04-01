/**
 * @file    version.h
 * @brief   Library version.
 */

#ifndef CRY_VERSION_H_
#define CRY_VERSION_H_

/** Library version major number */
#define CRY_MAJOR   0
/** Library version minor number */
#define CRY_MINOR   1
/** Library version patch number */
#define CRY_PATCH   1
/** Library version number */
#define CRY_VERSION (((CRY_MAJOR << 8) | (CRY_MINOR << 8)) | CRY_PATCH)

/**
 * Compile-time library version.
 *
 * For consistency, should be equal to the CRY_VERSION value.
 *
 * @return    compile time library version.
 */
long cry_version(void);

#endif /* CRY_VERSION_H_ */
