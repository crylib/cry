CC := gcc

CFLAGS := -O0 -g3 -Wall
CPPFLAGS = -I. -I../include -I../src
LDFLAGS := --coverage -Wl,--wrap=malloc,--wrap=realloc,--wrap=calloc,--wrap=free

libcry := ../libcry.a

objects := \
	test.o \
	malloc_mock.o \
	version_test.o \
	memxor_test.o \
	base64_test.o \
	mpi_test.o \
	des_test.o \
	aes_test.o \
	trivium_test.o \
	arc4_test.o \
	hill_test.o \
	affine_test.o \
	crc_test.o \
	hmac_test.o \
	cmac_test.o \
	rsa_test.o \
	md5_test.o \
	sha1_test.o \
	sha256_test.o \
	sha512_test.o \
	ecp_test.o \
	elgamal_test.o \
	dsa_test.o \
	ecdsa_test.o \
	utils_test.o

#	rand_test.o \
	dh_test.o \
	ecdh_test.o \

.PHONY: all clean 

all: $(objects) $(libcry)
	$(CC) $(LDFLAGS) -o test $^

clean:
	$(RM) $(objects) test
	$(RM) `find . -type f \( -name \*.gcda -o -name \*.gcno -o -name \*.gcov \)`

$(libcry):
	$(MAKE) -C ..

$(objects): ../include/cry/config.h

coverage: test
	@mkdir -p out/coverage
	@lcov -q -z -d ../build
	@lcov -q -c -i -d ../build -o base.info
	@./test $(TESTS)
	@lcov -q -c -d ../build -o cry.info
	@lcov -q -a base.info -a cry.info -o cry.info
	@genhtml -q -o out/coverage cry.info 
	@rm base.info cry.info

valgrind: all
	@mkdir -p out/valgrind
	valgrind --log-file=val.log --leak-check=full --show-leak-kinds=all --track-origins=yes ./test $(TESTS)
	@echo "<html><pre>" > valgrind.html
	@cat val.log >> valgrind.html
	@echo "</pre></html>" >> valgrind.html
	@mv valgrind.html out/valgrind/index.html
	@rm val.log

