/*
 * Copyright (c) 2013-2017, Davide Galassi. All rights reserved.
 *
 * This file is part of CRY software.
 *
 * CRY is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with CRY; if not, see <http://www.gnu/licenses/>.
 */

#include "mpi_pvt.h"

/*
 * Generate a random odd number of the desired bit length from a
 * pseudo-random source. Generally the output of the rand number generator
 * will have the first and the last bits set. Setting the last big ensures
 * the number is odd; no even numbers are primes. Setting the first bit
 * ensures that the generated number really is of the desired bit length.
 *
 * When generating RSA keys, people usually set the first two bits of all
 * their potential primes. That way, if you multiply two primes of the
 * same bit length together, they'll produce a result that's exactly twice
 * the bit length. When people talk about the "bit length of an RSA key",
 * they're generally talking about the size of such a product.
 *
 * To determine whether a number is prime, we use the Rabin-Miller test,
 * which can determine primality with high probability. Every time you run
 * the Rabin-Miller test and the test reports the number "may be prime,"
 * the actual probability of the number being prime increases dramatically.
 * By the time you've run five iterations and have received "may be prime"
 * every time, the odds of the random value's not being prime aren't worth
 * worrying about.
 *
 * If you are generating a prime number for use in Diffie-Hellman key exchange
 * (i.e., a "safe" prime), you should test the extra conditions before you
 * even check to see if the number itself is prime because doing so will
 * speed up tests.
 *
 * In this code, we explicitly attempt division for the first 100 primes,
 * although we recommend trying more primes than that (OpenSSL itself tries
 * 2,048, a widely recommended number.)
 */

static const int small_primes[] =
{
      3,    5,    7,   11,   13,   17,   19,   23,
     29,   31,   37,   41,   43,   47,   53,   59,
     61,   67,   71,   73,   79,   83,   89,   97,
    101,  103,  107,  109,  113,  127,  131,  137,
    139,  149,  151,  157,  163,  167,  173,  179,
    181,  191,  193,  197,  199,  211,  223,  227,
    229,  233,  239,  241,  251,  257,  263,  269,
    271,  277,  281,  283,  293,  307,  311,  313,
    317,  331,  337,  347,  349,  353,  359,  367,
    373,  379,  383,  389,  397,  401,  409,  419,
    421,  431,  433,  439,  443,  449,  457,  461,
    463,  467,  479,  487,  491,  499,  503,  509,
    521,  523,  541,  547,  557,  563,  569,  571,
    577,  587,  593,  599,  601,  607,  613,  617,
    619,  631,  641,  643,  647,  653,  659,  661,
    673,  677,  683,  691,  701,  709,  719,  727,
    733,  739,  743,  751,  757,  761,  769,  773,
    787,  797,  809,  811,  821,  823,  827,  829,
    839,  853,  857,  859,  863,  877,  881,  883,
    887,  907,  911,  919,  929,  937,  941,  947,
    953,  967,  971,  977,  983,  991,  997, -103
};


/*
 * Try simple division with all our small primes. This is, for each prime,
 * if it evenly divides p, return 0. Note that this obviously doesn't work
 * if we're checking a prime number that's in the list!
 *
 * On error returns -1
 */
static int is_obviously_not_prime(const cry_mpi *p)
{
    int i, res = 0;
    cry_mpi d, m;
    cry_mpi_digit dig;

    d.data = &dig;
    d.used = 1;
    d.alloc = 1;
    d.sign = 0;

    if ((res = cry_mpi_init(&m)) < 0)
        return res;

    for (i = 0;  small_primes[i] > 0;  i++) {
        dig = (cry_mpi_digit)small_primes[i];
        if ((res = cry_mpi_mod(&m, p, &d)) < 0)
            break;
        res = cry_mpi_is_zero(&m);
        if (res)
            break;
    }
    cry_mpi_clear(&m);
    return res; /* More tests needed */
}


/*
 * b is how many times does 2 divide p - 1. This gets returned.
 * m is (p-1)/(2^b).
 */
static int calc_b_and_m(cry_mpi *x, const cry_mpi *p)
{
    int ret = 0;
    cry_mpi one;
    cry_mpi_digit dig = 1;

    one.data = &dig;
    one.used = 1;
    one.alloc = 1;
    one.sign = 0;

    if (cry_mpi_init_copy(x, p) < 0)
        return -1;

    if ((ret = cry_mpi_sub(x, x, &one)) < 0)
        goto e;

    for (ret = 0; !cry_mpi_is_odd(x); ret++) {
        if (cry_mpi_shr(x, x, 1) < 0) { /* div by 2 */
            ret = -1;
            break;
        }
    }
e:  if (ret < 0)
        cry_mpi_clear(x);
    return ret;
}

/* 0 <= r < max */

static int passes_miller_rabin(const cry_mpi *p)
{
    int res;
    cry_mpi a, m, z, tmp;
    cry_mpi one;
    cry_mpi_digit dig = 1;
    unsigned int b, i;

    if ((res = cry_mpi_init_list(&a, &m, &z, &tmp, NULL)) < 0)
        return res;

    one.data = &dig;
    one.used = 1;
    one.alloc = 1;
    one.sign = 0;

    b = calc_b_and_m(&m, p);
    if ((int)b < 0) {
        res = (int)b;
        goto e;
    }

    /* Initialize a as a random number less than p */
    if ((res = cry_mpi_rand_range(&a, p)) < 0)
        goto e;

    /* z = a^m mod p */
    if ((res = cry_mpi_mod_exp(&z, &a, &m, p)) < 0)
        goto e;

    /* If z = 1 at the start, pass! */
    if (z.used == 1 && z.data[0] == 1) {
        res = 1;
        goto e;
    }

    for (i = 0; i < b; i++) {
        if (z.used == 1 && z.data[0] == 1) {
            res = 0;
            goto e;
        }

        /* if z = p - 1, pass! */
        if ((res = cry_mpi_copy(&tmp, &z)) < 0)
            goto e;
        if ((res = cry_mpi_add(&tmp, &tmp, &one)) < 0)
            goto e;
        if ((res = cry_mpi_cmp(&tmp, p)) == 0) {
            res = 1;
            goto e;
        }

        /* z = z^2 mod p */
        if ((res = cry_mpi_mul(&z, &z, &z)) < 0)
            goto e;
        if ((res = cry_mpi_mod(&z, &z, p)) < 0)
            goto e;
    }

    /* If z = p-1, pass! */
    if ((res = cry_mpi_add(&z, &z, &one) < 0))
        goto e;
    res = (cry_mpi_cmp(&z, p) == 0) ? 1 : 0;
e:  cry_mpi_clear_list(&a, &m, &z, &tmp, NULL);
    return res;
}

#define MILLER_ITER_NO  5

int cry_mpi_is_prime(const cry_mpi *p)
{
    int i, res;

    if ((res = is_obviously_not_prime(p)) < 0)
        return res;
    else if (res == 1)
        return 0;
    /* Is not obviously prime, proceed to Miller Rabin test */
    for (i = 0; i < MILLER_ITER_NO; i++) {
        if ((res = passes_miller_rabin(p)) <= 0)
            break; /* -1 on error; 0 if test is not passed */
    }
    return res;
}

/* Avoid infinite loops in bad machines */
#define ITERMAX  5000

/*
 * Assumes we only ever want to generate primes with the number of bits
 * multiple of 8.
 */
int cry_mpi_prime(cry_mpi *p, unsigned int bits, unsigned int *iter)
{
    int res = -1;
    unsigned int i, itermax;

    if (bits & 0x07 || bits == 0)
        return -1; /* Not a multiple of 8 bit */

    itermax = (iter) ? *iter : ITERMAX;
    i = 0;
    while (i < itermax) {
        i++;
        if ((res = cry_mpi_rand(p, bits)) < 0)
            break;
        /*
         * Set the least significant bit and the two most significant bits.
         * Lsb is set because any prime is odd.
         * Msb is set to ensure that the prime has the desired length.
         * The msb - 1 is set because that way when you multiply two big
         * primes of the same bit length together they'll produce a result
         * that is exactly twice their bit length.
         * This last property is especially usefull for RSA application.
         */
        cry_mpi_set_bit(p, 0);
        cry_mpi_set_bit(p, bits - 1);
        cry_mpi_set_bit(p, bits - 2);

        if ((res = cry_mpi_is_prime(p)) < 0) {
            break;
        } else if (res == 1) {
            res = 0;
            break;
        }
    }
    if (i == itermax)
        res = -1;
    if (iter)
        *iter = i;
    return res;
}

