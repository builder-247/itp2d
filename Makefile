# Add some standards-compliance and warning flags to the compiler
flags := $(CXXFLAGS) -Isrc -Wall -Wextra -std=c++98 -pedantic \
	-Wconversion -Wshadow -Wpointer-arith

# Add the version information if it is available
version := $(shell git describe --always --dirty)
ifdef version
	flags += -DVERSION=\"$(version)\"
endif

# Choose linear algebra library. Default linear algebra library is just plain
# CBLAS & LAPACK
# Possible values: plain, MKL, ACML
lalib := plain

# Default command for generating HTML documentation from Markdown
markdown := markdown2

# Override with local modifications written in Makefile local.mk, if such file
# exists
-include local.mk

# Add debug information with the binary. You can remove this if you want
# smaller binaries
flags += -g

# Uncomment this if you want to squeeze some more speed, although the
# difference is minimal
#flags += -DNDEBUG -O3 -ffast-math

# gtest and hdf5 are not that pedantic about standards so we need to silence a
# few warnings
flags += -Wno-long-long -Wno-variadic-macros

# Flags for FFTW and HDF5. Adjust if you have them installed somewhere away from the 
# compiler's search path
fftw_flags := -lfftw3
# Example: fftw_flags := -fftw3 -I/dir/where/fftw/include -L/dir/where/fftw/lib
hdf5_flags := -lhdf5 -lhdf5_cpp

# Flags for the linear algebra library
ifeq ($(lalib), MKL)
blas_flags := -lmkl
lapack_flags := -lmkl -lmkl_lapack64
flags += -DUSE_MKL
else
ifeq ($(lalib), ACML)
blas_flags := -lcblas -lacml
lapack_flags := -lacml
else
ifeq ($(lalib), plain)
blas_flags := -lcblas
lapack_flags := -llapack -lblas
else
$(error Variable lalib has invalid value "$(lalib)")
endif
endif
endif

# All libraries
lib_flags := -fopenmp $(fftw_flags) $(hdf5_flags) $(blas_flags) $(lapack_flags)

# Path to unpacked gtest source
gtest_dir = src/gtest
gtest_headers = $(gtest_dir)/include/gtest/*.h \
                $(gtest_dir)/include/gtest/internal/*.h
gtest_srcs = $(gtest_dir)/src/*.cc $(gtest_dir)/src/*.h $(gtest_headers)

# Flags for gtest
test_flags := -I $(gtest_dir)/include -I $(gtest_dir)
test_lib_flags := -pthread

progs := itp2d run_tests
src := $(wildcard src/*.cpp)
hdr := $(wildcard src/*.hpp)
obj := $(patsubst src/%.cpp,obj/%.o,$(src)) obj/gtest-all.o
dep := $(patsubst obj/%.o,.deps/%.o.d,$(obj)) .deps/itp2d.d .deps/run_tests.d

# Make targets and rules follow

.PHONY: default check doc all clean depend internalchecks

default: itp2d

check: run_tests
	./run_tests

doc: README.html

all: itp2d run_tests README.html

clean:
	rm -f $(progs) $(obj) $(dep)

depend: $(dep)

# Check that git is available and that the working tree is not dirty
internalchecks:
ifndef version
	$(warning Warning: could not query git for the program version)
	$(info If the version could be determined, it would be saved to the binary \
		and all created datafiles in order to aid debugging)
endif
ifneq (,$(findstring dirty,$(version)))
	$(warning Warning: According to git, you are compiling from a dirty working directory)
endif

data:
	mkdir -p data

.deps:
	mkdir -p .deps

obj:
	mkdir -p obj

.deps/%.o.d: src/%.cpp scripts/depmunger.py | .deps
	$(CXX) -MM -MT obj/$*.o $< -MF $@

.deps/itp2d.d: src/itp2d.cpp scripts/depmunger.py | .deps
	$(CXX) -MM -MT itp2d $< | scripts/depmunger.py > $@

.deps/run_tests.d: src/run_tests.cpp scripts/depmunger.py | .deps
	$(CXX) -MM -MT run_tests $< | scripts/depmunger.py > $@

obj/test_%.o: src/test_%.cpp $(gtest_headers)| obj
	$(CXX) $(flags) $(test_flags) -fopenmp -c $< -o $@

obj/run_tests.o: src/run_tests.cpp $(gtest_headers) | obj
	$(CXX) $(flags) $(test_flags) -fopenmp -c $< -o $@

obj/gtest-all.o: $(gtest_srcs) | $(gtest_dir)
	$(CXX) $(CXXFLAGS) $(test_flags) -c $(gtest_dir)/src/gtest-all.cc -o $@

obj/%.o: src/%.cpp | obj
	$(CXX) $(flags) -fopenmp -c $< -o $@

itp2d: obj/itp2d.o .deps/itp2d.d | data internalchecks
	$(CXX) $(flags) $(filter %.o,$^) $(lib_flags) -o $@

run_tests: obj/run_tests.o obj/gtest-all.o .deps/run_tests.d | data internalchecks
	$(CXX) $(flags) $(filter %.o,$^) $(lib_flags) $(test_lib_flags) -o $@

%.html: %.mdwn
	$(markdown) $< > $@

# include dependencies

-include $(dep)
