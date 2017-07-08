# src/subdir.mk

objs-y := version.o \
          memxor.o \
          base64.o \
          des.o \
          aes.o \
		  aes_wrap.o \
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
		  ecdsa.o

subdirs-y += mpi crc prng sum ecp

