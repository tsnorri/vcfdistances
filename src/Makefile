include ../local.mk
include ../common.mk

.PRECIOUS: vcf_reader.cc

OBJECTS		=	calculate_distances.o \
				cmdline.o \
				main.o

all: vcfdistances

clean:
	$(RM) $(OBJECTS) vcfdistances cmdline.c cmdline.h

vcfdistances: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS)

main.c: cmdline.c
cmdline.c: config.h

# Use the name config.h since cmdline.c includes it automatically.
config.h: SHELL := /bin/bash
config.h: Makefile ../.git
	IFS=$$'\n'; \
	TAGS=`git tag -l --points-at HEAD`; \
	HASH=`git rev-parse --short --verify HEAD`; \
	VERSION="DEV"; \
	if [ -n "$${TAGS}" ]; \
	then \
		for name in $${TAGS}; \
		do \
			if [[ "$${name}" == v* ]]; \
			then \
				VERSION="$${name:1}"; \
				break; \
			fi; \
		done; \
	fi; \
	printf "#define CMDLINE_PARSER_VERSION \"%s %s\"\n" "$${VERSION}" "$${HASH}" > config.h
