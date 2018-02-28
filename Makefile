include local.mk
include common.mk

DEPENDENCIES = lib/libbio/src/libbio.a
ifeq ($(shell uname -s),Linux)
	DEPENDENCIES	+= 
endif


.PHONY: all clean-all clean clean-dependencies dependencies

all: dependencies
	$(MAKE) -C src

clean-all: clean clean-dependencies

clean:
	$(MAKE) -C src clean

clean-dependencies:
	$(RM) -r lib/swift-corelibs-libdispatch/libdispatch-build

dependencies: $(DEPENDENCIES)
	
lib/libbio/src/libbio.a:
	cd lib/libbio && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	BOOST_INCLUDE="$(BOOST_INCLUDE)" \
	$(MAKE)
	
lib/swift-corelibs-libdispatch/libdispatch-build/libdispatch.a:
	
