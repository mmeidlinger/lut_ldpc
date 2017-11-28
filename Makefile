# This Makefile builds itpp and ldpc-lut


# ========== Set options for LDPC_LUT 
SRCDIR := src
PROGDIR := prog
DOCDIR := doc
BUILD_DIR := build
# Binaries and headers for external resources
LIBDIR := lib
INCDIR := include

LIB_PATH := -L$(LIBDIR)
# Set libraries
# Set libraries. Older gcc version do not autocomplete missing libraries by linking them dynamically, so if we compile on systens using older gcc versions, we need to specify the explicityl
LIBS_DYNAMIC :=   -llapack -lfftw3 -lblas -lgomp -lpthread -pthread  -static-libstdc++ 
LIBS_STATIC := -Wl,-Bstatic -lboost_filesystem -lboost_program_options -lboost_system -lglpk
LIBS_GCC4_DYNAMIC := -lamd -lcolamd -lz -lltdl -ldl
LIBS := $(LIBS_DYNAMIC) $(LIBS_GCC4_DYNAMIC) $(LIBS_STATIC)

CFLAGS := -std=c++11 -Wall
CXX := g++
LFLAGS :=


INCLUDE := -isystem $(INCDIR) -Isrc
DEBUG := 
BUILD_TYPE ?= Debug

# include bin/ directory in the path
#export PATH := bin:$(PATH)




# ========== Set options for ITPP 
ifeq ($(BUILD_TYPE),Release)
	EXECDIR := $(BUILD_DIR)/release/programs
	OBJDIR := $(BUILD_DIR)/release/objects
	ITPP_LIBNAME := libitpp_static.a
	LITPP := -litpp_static
 	CFLAGS += -DNDEBUG -O3
endif
ifeq ($(BUILD_TYPE),Debug)
	EXECDIR := $(BUILD_DIR)/debug/programs
	OBJDIR := $(BUILD_DIR)/debug/objects
	DEBUG += -ggdb
	ITPP_LIBNAME := libitpp_static_debug.a
	LITPP := -litpp_static_debug
endif

ITPP_TARGET= $(LIBDIR)/$(ITPP_LIBNAME)
ITPP_CMAKE_ARGS=-DITPP_SHARED_LIB=off -DHTML_DOCS=off -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DCMAKE_INSTALL_PREFIX=$(shell pwd)
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


NODEPS := clean clean-progs 

# ========== Degugging of Makefile
#$(info $$TARGETS is [${TARGETS}])
#$(info $$OSOURCES is [${OSOURCES}])
#$(info $$TSOURCES is [${TSOURCES}])

#================ Makefile targets

# Phony Targets
.PHONY: clean all doc

#TARGETS =  $(EXECDIR)/ber_sim $(EXECDIR)/de_sim $(EXECDIR)/degree_opt $(EXECDIR)/dat2alist $(EXECDIR)/ens2deg $(EXECDIR)/gen_code $(EXECDIR)/gen_ensemble $(EXECDIR)/alist2fullrank


# Compile
all: $(TARGETS)

# Targets: 
# Test explicit linking, for reference
#$(EXECDIR)/ber_sim: $(OBJECTS)
#	@echo " Linking $@..."
#	@mkdir -p $(dir $@)
#	$(CXX)  -o $@ \
#		$(OBJDIR)/$(notdir $@).o  \
#		$(OBJDIR)/LDPC_DE.o \
#		$(OBJDIR)/LDPC_Code_LUT.o \
#		$(OBJDIR)/LDPC_BER_Sim.o \
#		$(OBJDIR)/gitversion.o \
#		$(LIB_PATH) \
#		$(LITPP) \
#		-llapack \
#		-lfftw3 \
#		-lblas \
#		-lgomp \
#		-static-libstdc++ \
#		-Wl,-Bstatic \
#			-lboost_filesystem \
#			-lboost_program_options \
#			-lboost_system \
		

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
$(OBJDIR)/%.d: $(SRCDIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) $(INCLUDE) -MM $< > $@
	
# Create dependencies for cpp files with main function
$(OBJDIR)/%.d: $(PROGDIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) $(INCLUDE) -MM $< > $@
	
	
# Write Git Version to file for reproducibility of results
$(SRCDIR)/gitversion.cpp: ../.git/HEAD ../.git/index
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
	
# Installation	
install-debug:
	mkdir -p bin
	cp $(BUILD_DIR)/debug/programs/* bin/
	
install-release:
	mkdir -p bin
	cp $(BUILD_DIR)/release/programs/* bin/

install-remote:
	rsync -avzth --delete $(BUILD_DIR)/release/programs/ gate:~/epfl-tuwien-ldpc/cpp/bin/ 

# Delete everything 
clean:
	rm -rf $(BUILD_DIR)
	@rm -rf $(shell find $(DOCDIR)/* -not -name .gitignore -not -name Doxyfile -not -name README.md );
	rm -rf lib/libitpp*
	rm -rf $(INCDIR)/itpp	
	rm -rf $(SRCDIR)/gitversion.cpp
	
# Delete only executables (force relinking)
clean-progs:
	rm -rf $(EXECDIR);


