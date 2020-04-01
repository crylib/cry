[![Build Status](https://travis-ci.com/crylib/cry.svg?branch=master)](https://travis-ci.com/crylib/cry)
[![Coverage Status](https://coveralls.io/repos/github/crylib/cry/badge.svg?branch=master)](https://coveralls.io/github/crylib/cry?branch=master)


CRY
===

CRY is a small, comprehensive and portable cryptographic library prividing
a good collection of ciphers, hash functions, pseudo-random number generators
and multi-precision integers algorithms.

API Doxygen documentation [here](https://crylib.gitlab.io/cry/).


Symmetric ciphers
-----------------

### Block ciphers

- AES (Rijndael)
- DES and Triple DES

### Block ciphers mode of operation

- ECB (electronic codebook)
- CBC (cipher block chaining)
- CFB (cipher feedback)
- CTR (counter)
- GCM (Galois counter mode)

### Stream ciphers

- ARC4
- Trivium


Public key algorithms
---------------------

### Ciphers

- RSA (PKCS1 v1.5 padding)

### Secret exchange

- Diffie-Hellman (DH)
- ECDH

### Digital signature

- RSA (PKCS1 v2.1 PSS padding)
- DSA
- ECDSA
- Elgamal


Elliptic Curve
--------------

- EC group arithmetic.
- NIST recommended elliptic curve domain parameters over Fp.
- Brainpool standard curve domain parameters over Fp (RFC 5639).


Pseudo random numbers generators
--------------------------------

- AES-CTR CSPRNG
- LFSR-113 fast PRNG

Message authentication code
---------------------------

- HMAC
- CMAC


Hash
----

### Secure hash

- MD5
- SHA-1
- SHA-256
- SHA-512

### CRC

- CRC16-CCITT
- CRC16-IBM
- CRC16-DNP
- CRC32-Ethernet


Classical ciphers
-----------------

### Substitution ciphers

- Hill cipher
- Polyalphabetic Affine cipher


Multiple precision integers
---------------------------

- Basic arithmentic (add,sub,mul,div,mod,abs,exp,sqr)
- Modular exponentiation
- Modular inverse (Euclidean)
- GCD and LCM
- Probabilistic prime numbers generator (Miller-Rabin)
- Random mpi generator
- Comba, Karatsuba and Toom-Cook-3 multipliers


Utilities
---------

- AES Key Wrap (RFC-3394)
- Base64 encoder/decoder
- Memxor
