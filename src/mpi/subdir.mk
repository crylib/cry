# src/mpi/subdir.mk

objects-y := \
    mpi_core.o \
    mpi_list.o \
    mpi_abs.o \
    mpi_cmp.o \
    mpi_cmp_abs.o \
    mpi_add.o \
    mpi_add_abs.o \
    mpi_sub.o \
    mpi_sub_abs.o \
    mpi_mul.o \
    mpi_mul_abs.o \
    mpi_mul_baseline.o \
    mpi_mul_comba.o \
    mpi_mul_karatsuba.o \
    mpi_mul_toom3.o \
    mpi_div.o \
    mpi_div_abs.o \
    mpi_sqr.o \
    mpi_sqr_baseline.o \
    mpi_sqr_comba.o \
    mpi_sqr_karatsuba.o \
    mpi_sqrt.o \
    mpi_shl.o \
    mpi_shr.o \
    mpi_bin.o \
    mpi_str.o \
    mpi_exp.o \
    mpi_mod_exp.o \
    mpi_gcd.o \
    mpi_lcm.o \
    mpi_inv.o \
    mpi_rand.o \
    mpi_prime.o \
    mpi_print.o