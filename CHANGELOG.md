Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com).

Type of changes

* Added: for new features.
* Changed: for changes in existing functionality.
* Deprecated: for soon-to-be removed features.
* Removed: for now removed features.
* Fixed: for any bug fixes.
* Security: in case of vulnerabilities.

This project adheres to [Semantic Versioning](http://semver.org).

Given a version number MAJOR.MINOR.PATCH
* MAJOR incremented for incompatible API changes
* MINOR incremented for new functionalities
* PATCH incremented for bug fixes


[0.1.1] - Unreleased
--------------------

Added

* ARC4 stream cipher
* Elgamal digital signature

Fixed

* MPI "is_prime" function was not working for primes < 1010 (decimal)


[0.1.0] - 2019-09-23
--------------------

Added

* RSA sign/verify functions
* Schoolbook RSA
* AES Key Wrap (RFC-3394)

Changed

* RSA padding schemes correction
* MPI random prime generation sequence independence from MPI digit size
* Weak PRNG properly renamed to LFSR-113

Fixed

* MPI bit-set macro was setting the wrong bit when using 64-bit digits


[0.0.10] - 2019-06-08
---------------------

Added

* Some NIST recommended elliptic curve domain parameters over Fp.
* Some Brainpool standard curve domain parameters over Fp (RFC 5639).
* Malloc mock and Valgrind used to execute CI unit-tests>
* Uncrustify and clang-format project files.
* Stricter compiler warnings policy wrt size types and int sign conversions.


[0.0.9] - 2019-04-22
--------------------

Added

* SHA1 message digest
* Continuous code quality inspection.
  Tools: cppcheck, vera++, valgrind, gcovr, sonar-scanner, clang.

Changed

* Configurable MPI digit size and memory alloc quantum.
* MPI exponentiation refactory.

Fixed

* MPI product of negative value by zero returning negative zero.
* MPI potential uninitialized memory access when digit size is 1.
* MPI square unhandled integer overflow when digit size is 1.


[0.0.8] - 2019-03-16
--------------------

Added

* CFB block cipher mode of operation.
* CRC16-DNP.
* RSA keypair generator.

Changed

* Incremental CMAC computation.

Removed

* Trivial LRC and 8-bit checksum.
* System dependent PRNG.

Fixed

* MPI and RSA not initialized values.


[0.0.7] - 2018-12-26
--------------------

Added

* MPI Toom-Cook-3 multiplier.
* Trivium stream cipher.
* Hill cipher.
* Polyalphabetic affine cipher.

Changed

* Project released under the MIT license.


[0.0.6] - 2017-07-09
--------------------

Added

* ECP core arithmetic (add,dbl,mul).
* Load NIST-P256 EC parameters.
* DSA end ECDSA digital signature.
* ECP multiplication with sliding window algorithm.
* ECDH secret exchange.


[0.0.5] - 2016-10-02
--------------------

Added

* Diffie-Hellman secret exchange.
* OS agnostic CSPRNG with AES-CTR.
* CRC32 Ethernet.


[0.0.4] - 2016-05-07
--------------------

Added

* MPI Comba multiplier.
* MPI Karatsuba multiplier.
* MPI fast division.
* MPI square algorithm.
* HMAC with generic hash function.


[0.0.3] - 2015-11-28
--------------------

Added

* MPI abs, shl/shr, gcd, lcm, exp, mod-exp.
* MPI load/store from string and binary buffer.
* MPI inverse (Euclidean algorithm).
* Unit test framework.
* Improved build system using non recursive make.
* CSPRNG wrapper around OS specific PRNG..
* RSA (PKCS#1 v1.5).


[0.0.2] - 2014-11-22
--------------------

Added

* LRC.
* 8-bit trivial checksum.
* CRC16-IBM.
* MPI add, sub, mul, div.
* MPI print.


[0.0.1] - 2014-05-10
--------------------

Added

* AES block cipher.
* DES block cipher.
* CBC block cipher mode of operation.
* GCM block cipher mode of operation.
* CTR block cipher mode of operation.
* MD5 message digest.
* SHA256 message digest.
* CMAC message authentication code.
* Base64 encoder.
* CRC16-CCITT.
