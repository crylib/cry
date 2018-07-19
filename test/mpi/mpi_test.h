#ifndef TEST_MPI_MPI_TEST_H_
#define TEST_MPI_MPI_TEST_H_

#include "test.h"
#include <cry/mpi.h>


#define MPI_BUF_LEN   4
extern cry_mpi *g_mpi_buf;
#define g_mpi0 (&g_mpi_buf[0])
#define g_mpi1 (&g_mpi_buf[1])
#define g_mpi2 (&g_mpi_buf[2])
#define g_mpi3 (&g_mpi_buf[3])

#define A32_VAL   ((long)0xc92e20b9)

extern const char g_a6400_bin[800];
extern const char g_b6400_bin[800];

extern const char g_a8_b8_add_bin[1];
extern const char g_a8_b256_add_bin[32];
extern const char g_a1024_a512_add_bin[128];
extern const char g_a1024_a1024_add_bin[129];
extern const char g_a264_a256_add_sign_bin[33];

extern const char g_a8_b8_sub_bin[1];
extern const char g_a512_a256_sub_bin[64];
extern const char g_a1024_a256_sub_bin[128];
extern const char g_a264_a256_sub_sign_bin[33];


extern const char g_a6400_b6400_mul_bin[1600];


void mpi_setup(void);
void mpi_teardown(void);
#define MPI_RUN(test) RUNX(test, mpi_setup, mpi_teardown)

/*
 * Subtests
 */
void mpi_core_test(void);
void mpi_cmp_test(void);
void mpi_abs_test(void);
void mpi_add_test(void);
void mpi_sub_test(void);
void mpi_shl_test(void);
void mpi_shr_test(void);
void mpi_mul_test(void);

#endif /* TEST_MPI_MPI_TEST_H_ */