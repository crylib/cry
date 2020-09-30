CONFIG := config/config.mk

include $(CONFIG)

CC := gcc
AR := ar
AWK := awk
CP := cp
RM := rm -rf

source_dir := src
config_mk := $(CONFIG)
config_h := include/cry/config.h

# $(call normalstr,str)
build_name = $(shell $(CC) -dumpmachine)
# convert to lowercase
binary_dir = build/$(build_name)

# call $(src_to_bin_dir,list)
src_to_bin_dir = $(patsubst $(source_dir)%,$(binary_dir)%,$1)

target = $(binary_dir)/libcry.a

.SUFFIXES:

warnings := -Wall -Wextra -Wshadow -Wstrict-prototypes -Wmissing-prototypes

#-Wpedantic
#-Wconversion
includes-y := -Iinclude -Isrc

cflags-y := -MMD -MP $(warnings)
cflags-$(CRY_COVERAGE) += --coverage
cflags-$(CRY_SHORT_ENUMS) += -fshort-enums
cflags-$(CRY_OMIT_FRAME_POINTER) += -fomit-frame-pointer
cflags-$(CRY_NO_STACK_PROTECTOR) += -fno-stack-protector

ifeq ($(CRY_DEBUG),y)
cflags-y += -O0 -g3
else
cflags-y += -DNDEBUG
ifeq ($(CRY_SMALL_SIZE),y)
cflags-y += -Os
else
cflags-y += -O3
endif
endif

objects-y :=
paths-y	:=
objects_list :=

define include_subdir
$(shell mkdir -p $(call src_to_bin_dir,$1))
subdirs-y :=
current := $1
objects-y :=
include $1/subdir.mk
objects_list += $$(if $$(objects-y),$$(addprefix $1/,$$(objects-y)))
paths-y += $$(current)
subdirs-y := $$(addprefix $$(current)/, $$(subdirs-y))
$$(foreach subdir, $$(subdirs-y),$$(eval $$(call include_subdir,$$(subdir))))
endef

$(eval $(call include_subdir,src))

objects = $(call src_to_bin_dir,$(objects_list))
depends = $(patsubst %.o,%.d,$(objects))

CPPFLAGS := $(includes-y)
CFLAGS   := $(cflags-y)
AFLAGS   := $(aflags-y)
LDFLAGS  := $(lflags-y)

################################################################################
# Targets
################################################################################

.PHONY: all clean config test testclean doc

# Force serial run
all:
	$(MAKE) config
	$(MAKE) $(target)

clean:
	@echo "Cleanup ..."
	@$(RM) $(binary_dir) *.a
	@$(RM) `find . -type f \( -name \*.gcda -o -name \*.gcno \)`

config: $(config_mk)
	@printf "/*\n * Automatically generated from \"$^\".\n */\n\n" > tmp
	@$(AWK) -F= 'NF > 1 && $$1 !~ /^[# ]/ { print "#define", $$1; }' < $^ >> tmp
	@cmp -s tmp $(config_h) || (echo "Configuration update"; cp tmp $(config_h))
	@$(RM) tmp
	@echo ">>> Config : $(config_mk)"
	@cat $(config_h) | grep CRY_ | $(AWK) '{ printf(" * %s\n", $$2); }'

$(target): $(objects)
	$(AR) rcs $@ $^
	$(CP) $(target) .

$(objects): Makefile $(config_h)

$(config_h): $(config_mk)
	touch $(config_h)

$(binary_dir)/%.o: $(source_dir)/%.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

test: all
	$(MAKE) -C test

testclean:
	$(MAKE) -C test clean

doc:
	cd doc/doxy; ./build.sh

# include the .d files, if they do not exists then are generated by
# the rule the pattern rule.
# The `-` prefix suppress the warnings relative to the fact that initially
# make do not finds the files (before trying to create them).
ifneq ($(MAKECMDGOALS),clean)
-include $(depends)
endif

################################################################################
# Code quality tools reports targets
################################################################################

# Overwrite from command line with: `make SONAR_SCANNER=<path>`
SONAR_SCANNER := sonar-scanner

# Vera style checker rules to skip
VERA_SKIP := "T0(10|11|12|13|19)"

# CPPcheck: a general purpose static code checker
cppcheck:
	@echo "Running cppcheck ..."
	@cppcheck --language=c --std=c89 --force --enable=all --suppress=variableScope --suppress=unusedFunction --suppress=missingIncludeSystem --xml -I src -I include src 2> build/cppcheck-report.xml

# Vera++: static code checker focusing on code style issues
vera:
	@echo "Running vera++ ..."
	@find src include -type f -regextype sed -iregex ".*/*\.\(c\|h\)" -print | vera++ - -showrules -nodup 2> vera.tmp
	@cat vera.tmp | grep -v -E $(VERA_SKIP) | scripts/vera2report.perl > build/vera-report.xml
	@rm vera.tmp

# Valgrind: memory leaks and dynamic issues report
valgrind: test
	@echo "Running valgrind ..."
	@cd test; valgrind --xml=yes --xml-file=valgrind-report.xml --leak-check=full --show-leak-kinds=all --track-origins=yes ./test -v; cd ..
	@mv test/valgrind-report.xml build/valgrind-report.xml

# Gcovr: coverage report using unit tests (uses valgrind run)
coverage: valgrind
	@echo "Running gcovr ..."
	@gcovr --xml --root . > build/gcovr-report.xml

# Sonarqube: continuous inspection of code quality platform.
# Parse collected data and feed it into Sonarqube server.
sonar: cppcheck vera coverage valgrind
	@echo "Running sonar-scanner ..."
	@$(SONAR_SCANNER) -Dproject.settings=scripts/sonar-project.properties
