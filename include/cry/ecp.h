/**
 * @file    ecp.h
 * @brief   Elliptic curve point arithmetic
 */

#ifndef CRY_ECP_H_
#define CRY_ECP_H_

#include <cry/mpi.h>

#define CRY_ECP_GRP_SECP192R1   0
#define CRY_ECP_GRP_SECP224R1   1
#define CRY_ECP_GRP_SECP256R1   2
#define CRY_ECP_GRP_SECP384R1   3
#define CRY_ECP_GRP_SECP521R1   4
#define CRY_ECP_GRP_SECP192K1   5
#define CRY_ECP_GRP_SECP224K1   6
#define CRY_ECP_GRP_SECP256K1   7
#define CRY_ECP_GRP_BP256R1     8
#define CRY_ECP_GRP_BP384R1     9
#define CRY_ECP_GRP_BP512R1     10

/**
 * EC point
 */
struct cry_ecp {
    cry_mpi x;
    cry_mpi y;
    cry_mpi z;
};

typedef struct cry_ecp cry_ecp;

/*
 * EC group equation parameters.
 * y^2 = (x^3 + ax + b) % p
 */
struct cry_ecp_grp {
    cry_mpi p;  /**< Modulus */
    cry_mpi a;  /**< x coefficient */
    cry_mpi b;  /**< Constant */
    cry_ecp g;  /**< Generator point */
    cry_mpi n;  /**< Generator order */
};

typedef struct cry_ecp_grp cry_ecp_grp;


#ifdef __cplusplus
extern "C" {
#endif

int cry_ecp_init(cry_ecp *p);

void cry_ecp_clear(cry_ecp *p);

int cry_ecp_init_int(cry_ecp *p, long x, long y);

int cry_ecp_copy(cry_ecp *pr, const cry_ecp *p1);

int cry_ecp_add(cry_ecp *pr, const cry_ecp *p1, const cry_ecp *p2,
                const cry_ecp_grp *grp);

int cry_ecp_dbl(cry_ecp *pr, const cry_ecp *p1, const cry_ecp_grp *grp);

int cry_ecp_mul(cry_ecp *pr, const cry_ecp *p1, const cry_mpi *k,
                const cry_ecp_grp *grp);

int cry_ecp_init_list(cry_ecp *p, ... /* (cry_ecp *)NULL */);

void cry_ecp_clear_list(cry_ecp *p, ... /* (cry_ecp *)NULL */);

#ifdef __cplusplus
}
#endif

#define cry_ecp_is_zero(p) \
    cry_mpi_is_zero(&(p)->z)

#define cry_ecp_set_zero(p) do { \
    cry_mpi_set_int(&(p)->x, 1); \
    cry_mpi_set_int(&(p)->y, 1); \
    cry_mpi_set_int(&(p)->z, 0); \
} while (0)

#define cry_ecp_swap(a, b) do {     \
    cry_mpi_swap(&(a)->x, &(b)->x); \
    cry_mpi_swap(&(a)->y, &(b)->y); \
    cry_mpi_swap(&(a)->z, &(b)->z); \
} while (0)

int cry_ecp_grp_init(cry_ecp_grp *grp);

void cry_ecp_grp_clear(cry_ecp_grp *grp);

int cry_ecp_grp_load(cry_ecp_grp *grp, int grp_id);


#endif /* CRY_ECP_H_ */
