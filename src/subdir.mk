# src/subdir.mk

objects-y := \
    version.o \
    memxor.o \
    base64.o \
    des.o \
    aes.o \
    cbc.o \
    gcm.o \
    ctr.o \
    md5.o \
    sha256.o \
    cmac.o \
    hmac.o \
    rsa.o \
    dh.o \
    ecdh.o \
    dsa.o \
    ecdsa.o \
    trivium.o \
    misc.o

subdirs-y += mpi crc prng sum ecp classic