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

main.c : cmdline.c
