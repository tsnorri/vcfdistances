# Ignore Xcode's setting since the SDK may contain older versions of Clang and libc++.
unexport SDKROOT

WARNING_FLAGS	= -Wall -Werror -Wno-deprecated-declarations -Wno-unused -Wno-misleading-indentation
OPT_FLAGS		= -O2 -g

CFLAGS			= -std=c99   $(OPT_FLAGS) $(WARNING_FLAGS)
CXXFLAGS		= -std=c++1z $(OPT_FLAGS) $(WARNING_FLAGS)
CPPFLAGS		= -I../include -I../lib/libbio/include -I../lib/sdsl-lite/include $(BOOST_INCLUDE)
LDFLAGS			= ../lib/libbio/src/libbio.a $(LIBDISPATCH_LIBS) $(BOOST_LIBS) -lz

ifeq ($(shell uname -s),Linux)
	#CPPFLAGS    += -I../lib/swift-corelibs-libdispatch
	CPPFLAGS += -I../lib/libdispatch -I../lib/libpwq/include
endif

GENGETOPT		?= gengetopt
CMAKE			?= cmake
MKDIR			?= mkdir
NINJA			?= ninja
RAGEL			?= ragel


%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<"

%.cc: %.rl
	$(RAGEL) -L -C -G2 -o $@ $<

%.dot: %.rl
	$(RAGEL) -V -p -o $@ $<

%.pdf: %.dot
	$(DOT) -Tpdf $< > $@
