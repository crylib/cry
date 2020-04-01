################################################################################
# CRY build configuration.
#
# This file is included by Makefile.
################################################################################

################################################################################
# Compiler
################################################################################

CRY_DEBUG=y
#CRY_SMALL_SIZE=y
#CRY_OMIT_FRAME_POINTER=y
#CRY_NO_STACK_PROTECTOR=y
#CRY_COVERAGE=y

################################################################################
# Platform
################################################################################

#CRY_ARCH_X86_64=y
#CRY_ARCH_X86=y

################################################################################
# Algorithms tweaks
################################################################################

# Use Extended Euclidean Algorithm to compute inverse.
# If not defined the faster Binary Extended Algorithm (HAC 14.4.3) is used.
CRY_MPI_INV_EUCLID=y

# Small footprint division. Slower.
#CRY_MPI_DIV_SMALL=y

# Comba as baseline multiplier
CRY_MPI_MUL_COMBA=y

# Karatsuba multiplier enabled, see the mpi_mul.c to set the CUTOFF
CRY_MPI_MUL_KARATSUBA=y

# Toom-3 multiplier enabled, see the mpi_mul.c to set the CUTOFF
CRY_MPI_MUL_TOOM3=y

# Elliptic curve point multiplication with windowed algorithm
CRY_ECP_MUL_WIN=y

# Elliptic curve point multiplication with sliding window algorithm (faster)
# Have precedence against the simple windowed one.
CRY_ECP_MUL_SLIDING_WIN=y

################################################################################
# Elliptic curves group parameters
################################################################################

CRY_ECP_GRP_SECP192R1_ENABLE=y
CRY_ECP_GRP_SECP224R1_ENABLE=y
CRY_ECP_GRP_SECP256R1_ENABLE=y
CRY_ECP_GRP_SECP384R1_ENABLE=y
CRY_ECP_GRP_SECP521R1_ENABLE=y
CRY_ECP_GRP_SECP192K1_ENABLE=y
CRY_ECP_GRP_SECP224K1_ENABLE=y
CRY_ECP_GRP_SECP256K1_ENABLE=y
CRY_ECP_GRP_BP256R1_ENABLE=y
CRY_ECP_GRP_BP384R1_ENABLE=y
CRY_ECP_GRP_BP512R1_ENABLE=y

################################################################################
# Miscelanea configurations
################################################################################

#CRY_MPI_DEBUG_CONF=y
#CRY_ASSERT_CHECK=y

