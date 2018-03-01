include local.mk
include common.mk

DEPENDENCIES = lib/libbio/src/libbio.a
ifeq ($(shell uname -s),Linux)
	DEPENDENCIES    +=  lib/libdispatch/libdispatch-build/src/libdispatch.a
	#DEPENDENCIES	+=	lib/swift-corelibs-libdispatch/libdispatch-build/src/libdispatch.a
	DEPENDENCIES    +=  lib/libpwq/libpwq-build/libpthread_workqueue.a
endif


.PHONY: all clean-all clean clean-dependencies dependencies

all: dependencies
	$(MAKE) -C src

clean-all: clean clean-dependencies

clean:
	$(MAKE) -C src clean

clean-dependencies:
	$(RM) -r lib/swift-corelibs-libdispatch/libdispatch-build
	$(RM) -r lib/libdispatch/libdispatch-build
	$(RM) -r lib/libpwq/libpwq-build
	$(RM) -r lib/lemon/build/lemon/libemon.a

dependencies: $(DEPENDENCIES)
	
lib/libbio/src/libbio.a:
	cd lib/libbio && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	BOOST_INCLUDE="$(BOOST_INCLUDE)" \
	$(MAKE)
	
lib/swift-corelibs-libdispatch/libdispatch-build/src/libdispatch.a: 
	$(RM) -r lib/swift-corelibs-libdispatch/libdispatch-build && \
	cd lib/swift-corelibs-libdispatch && \
	$(MKDIR) libdispatch-build && \
	cd libdispatch-build && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	CFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CFLAGS)" \
	CXXFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CXXFLAGS)" \
	LDFLAGS="$(SYSTEM_LDFLAGS)" \
	$(CMAKE) -G Ninja -DBUILD_SHARED_LIBS=OFF -DWITH_BLOCKS_RUNTIME=/usr/lib ..
	$(NINJA) -C lib/swift-corelibs-libdispatch/libdispatch-build

lib/libdispatch/libdispatch-build/src/libdispatch.a: lib/libpwq/libpwq-build/libpthread_workqueue.a
	$(RM) -r lib/libdispatch/libdispatch-build && \
	cd lib/libdispatch && \
	$(MKDIR) libdispatch-build && \
	cd libdispatch-build && \
	CFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CFLAGS)" \
	CXXFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CXXFLAGS)" \
	LDFLAGS="$(SYSTEM_LDFLAGS)" \
	../configure --cc="$(CC)" --c++="$(CXX)" --release -- \
		-DPTHREAD_WORKQUEUE_INCLUDE_DIRS=../../libpwq/include \
		-DPTHREAD_WORKQUEUE_LIBRARIES=../../libpwq/libpwq-build/libpthread_workqueue.a \
		-DBLOCKS_RUNTIME_LIBRARIES=""
	$(MAKE) -C lib/libdispatch/libdispatch-build VERBOSE=1

lib/libpwq/libpwq-build/libpthread_workqueue.a:
	$(RM) -r lib/libpwq/libpwq-build && \
	cd lib/libpwq && \
	$(MKDIR) libpwq-build && \
	cd libpwq-build && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	CFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CFLAGS)" \
	CXXFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CXXFLAGS)" \
	LDFLAGS="$(SYSTEM_LDFLAGS)" \
	$(CMAKE) -DSTATIC_WORKQUEUE=ON ..
	$(MAKE) -C lib/libpwq/libpwq-build VERBOSE=1

