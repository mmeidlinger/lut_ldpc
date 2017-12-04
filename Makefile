# This Makefile builds itpp and ldpc-lut


# ========== Set options for LDPC_LUT 
SRCDIR := src
PROGDIR := prog
DOCDIR := doc
BUILD_DIR := build
INSTALLDIR := bin
# Binaries and headers for external resources
LIBDIR := lib
INCDIR := include
BUILD_TYPE ?= Release
LINK_TYPE ?= static

# ========== Set compiler based on Operating system
ifeq ($(OS),Windows_NT)
    echo "lut_ldpc is currently not supported for Windows. There is no fundamental reason it shouldn't work though, you just need to figure out the compilation and linking process yourself";
    exit
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
		CXX := g++
		LIBS_DYNAMIC :=   -llapack -lfftw3 -lblas -lgomp -lpthread -pthread  -lboost_filesystem -lboost_program_options -lboost_system 
		LIBS_STATIC :=  -llapack -lfftw3 -lblas -lgomp -lpthread -pthread  -static-libstdc++ -Wl,-Bstatic -lboost_filesystem -lboost_program_options -lboost_system
		LIB_DYN_EXTENSION := so
    endif
    ifeq ($(UNAME_S),Darwin)
        CXX := clang++
		LIBS_STATIC := -llapack -lfftw3 -lblas  -lboost_filesystem -lboost_program_options -lboost_system
		# OSX is not very well suited for static linking. Hence, boost is linked dynamically here
		LIBS_DYNAMIC :=  $(LIBS_STATIC)
		LIB_DYN_EXTENSION := dylib
    endif
endif

LIB_PATH := -L$(LIBDIR)
# Set libraries
ifeq ($(LINK_TYPE),static)
	LIBS := $(LIBS_STATIC)
endif
ifeq ($(LINK_TYPE),dynamic)
	LIBS := $(LIBS_DYNAMIC)
endif


CFLAGS := -std=c++11 -Wall
LFLAGS :=


INCLUDE := -isystem $(INCDIR) -Isrc
DEBUG := 





# ========== Set options for ITPP 
ifeq ($(BUILD_TYPE),Release)
	EXECDIR := $(BUILD_DIR)/release/programs
	OBJDIR := $(BUILD_DIR)/release/objects
 	CFLAGS += -DNDEBUG -O3
 	ifeq ($(LINK_TYPE),static)
		ITPP_LIBNAME := libitpp_static.a
		LITPP := -litpp_static
		LITPP_IS_SHARED := off
	endif
	ifeq ($(LINK_TYPE),dynamic)
		ITPP_LIBNAME := libitpp.$(LIB_DYN_EXTENSION)
		LITPP := -litpp
		LITPP_IS_SHARED := on
	endif
endif
ifeq ($(BUILD_TYPE),Debug)
	EXECDIR := $(BUILD_DIR)/debug/programs
	OBJDIR := $(BUILD_DIR)/debug/objects
	DEBUG += -ggdb
	ifeq ($(LINK_TYPE),static)
		ITPP_LIBNAME := libitpp_static_debug.a
		LITPP := -litpp_static_debug
		LITPP_IS_SHARED := off
	endif
	ifeq ($(LINK_TYPE),dynamic)
		ITPP_LIBNAME := libitpp_debug.$(LIB_DYN_EXTENSION)
		LITPP := -litpp_debug
		LITPP_IS_SHARED := on
	endif
endif

ITPP_TARGET= $(LIBDIR)/$(ITPP_LIBNAME)
ITPP_CMAKE_ARGS=-DITPP_SHARED_LIB=$(LITPP_IS_SHARED) -DHTML_DOCS=off -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DCMAKE_INSTALL_PREFIX=$(shell pwd)
ITPP_HEADERS := $(shell find itpp/itpp -name '*.h'  -o -name '*.hpp')
ITPP_SOURCES := $(shell find itpp/itpp -name '*.c' -o -name '*.cpp')




SOURCES_SRC := $(shell find $(SRCDIR) -type f -name *.cpp)
SOURCES_PROG := $(shell find $(PROGDIR) -type f -name *.cpp)
SOURCES = $(SOURCES_SRC) $(SOURCES_PROG)

OBJECTS_SRC := $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SOURCES_SRC:.cpp=.o))
OBJECTS_PROG := $(patsubst $(PROGDIR)/%,$(OBJDIR)/%,$(SOURCES_PROG:.cpp=.o))
OBJECTS :=  $(OBJECTS_SRC) $(OBJECTS_PROG)
OBJECTS += $(OBJDIR)/gitversion.o

TARGETS = $(patsubst $(OBJDIR)/%,$(EXECDIR)/%,$(OBJECTS_PROG:.o=))


DFILES := $(OBJECTS:.o=.d)


NODEPS := clean clean-all

# ========== Degugging of Makefile
#$(info $$TARGETS is [${TARGETS}])
#$(info $$OBJECTS is [${OBJECTS}])
#$(info $$OSOURCES is [${OSOURCES}])
#$(info $$TSOURCES is [${TSOURCES}])

#================ Makefile targets



# Phony Targets
.PHONY: clean clean-all all doc


# Compile
all: $(TARGETS)


# other target implicit linking
$(TARGETS): $(OBJECTS)
	@echo "Linking $@..."
	@mkdir -p $(dir $@)
	$(CXX)  -o $@ \
		$(OBJDIR)/$(notdir $@).o  \
		$(OBJECTS_SRC) \
		$(LIB_PATH) \
		$(LITPP) \
		$(LIBS)


ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
-include $(DFILES)          # slient fail include: if the *.d files are not existing, 
	  			      # make will not abort execution
endif

# Create dependencies for cpp files without main function
$(OBJDIR)/%.d: $(SRCDIR)/%.cpp  $(ITPP_TARGET)
	@mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) $(INCLUDE) -MM $< > $@
	
# Create dependencies for cpp files with main function
$(OBJDIR)/%.d: $(PROGDIR)/%.cpp  $(ITPP_TARGET)
	@mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) $(INCLUDE) -MM $< > $@
	
	
# Write Git Version to file for reproducibility of results
$(SRCDIR)/gitversion.cpp: .git/HEAD .git/index
	echo "const char *gitversion = \"$(shell git rev-parse HEAD)\";" > $@

# Build ITPP and install (headers go into include/, the lib itself into lib/)
$(ITPP_TARGET): $(ITPP_HEADERS) $(ITPP_SOURCES) 
	cd itpp && mkdir -p build && cd build && \
	cmake .. $(ITPP_CMAKE_ARGS) && \
	make && make install && cd ../.. ; \
	rm -r share && rm -r $(LIBDIR)/pkgconfig;

# Create object files without main function
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)/%.d $(ITPP_TARGET) 
	@mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) $(INCLUDE) $(DEBUG) -o $@ -c $<
	
# Create object files with main function
$(OBJDIR)/%.o: $(PROGDIR)/%.cpp $(OBJDIR)/%.d $(ITPP_TARGET) 
	@mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) $(INCLUDE) $(DEBUG)  -o $@ -c $<

# Documentation
doc:
	doxygen $(DOCDIR)/Doxyfile
	./doc/doxy_post.sh
	
# Installation	
install-debug:
	mkdir -p $(INSTALLDIR)
	cp $(BUILD_DIR)/debug/programs/* $(INSTALLDIR)/

	
install-release:
	mkdir -p $(INSTALLDIR)
	cp $(BUILD_DIR)/release/programs/* $(INSTALLDIR)/


install: install-release
	
# Delete everything 
clean-all:
	rm -rf $(BUILD_DIR);
	rm -rf doc/{html,latex}
	rm -rf lib/libitpp*
	rm -rf itpp/build*
	rm -rf $(INCDIR)/itpp	
	rm -rf $(SRCDIR)/gitversion.cpp
	
# Delete only executables (force relinking)
clean:
	rm -rf $(BUILD_DIR);


