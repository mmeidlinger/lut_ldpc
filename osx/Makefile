# Create dependencies. The actual build is done by xcode

ITPP_TARGET=lib/libitpp_static_debug.a
ITPP_CMAKE_ARGS=-DITPP_SHARED_LIB=off -DHTML_DOCS=off -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$(shell pwd)


ITPP_HEADERS := $(shell find ../itpp/itpp -name '*.h'  -o -name '*.hpp')
ITPP_SOURCES := $(shell find ../itpp/itpp -name '*.c' -o -name '*.cpp')

all: $(ITPP_TARGET) gitversion.cpp

# Build ITPP and install (headers go into include/, the lib itself into lib/)
$(ITPP_TARGET): $(ITPP_HEADERS) $(ITPP_SOURCES) 
	mkdir -p itpp_build_osx && cd itpp_build_osx && \
	cmake ../../itpp $(ITPP_CMAKE_ARGS) && \
	make && make install && cd .. ; \
	rm -rf share && rm -rf lib/pkgconfig;

# Write Git Version to file for reproducibility of results
gitversion.cpp: ../../.git/HEAD ../../.git/index
	echo "const char *gitversion = \"$(shell git rev-parse HEAD)\";" > $@

