# Add some standards-compliance and warning flags to the compiler
flags := $(CXXFLAGS) -Isrc -Wall -Wextra -std=c++98 -pedantic \
	-Wconversion -Wshadow -Wpointer-arith

# gtest and hdf5 are not that pedantic about standards so we need to silence a
# few warnings
flags += -Wno-long-long -Wno-variadic-macros

# Add flag for OpenMP
flags += -fopenmp

# Add the version information if it is available
version := $(shell git describe --always --dirty)
ifdef version
	flags += -DITP2D_VERSION=\"$(version)\"
endif

lib_flags := -fopenmp -lrt -lfftw3 -lhdf5 -lhdf5_cpp
inc_flags :=
# Path to unpacked gtest source
gtest_dir = src/gtest

# Override with local modifications written in Makefile local.mk, if such file
# exists
-include Makefile.local

gtest_headers = $(gtest_dir)/include/gtest/*.h \
                $(gtest_dir)/include/gtest/internal/*.h
gtest_srcs = $(gtest_dir)/src/*.cc $(gtest_dir)/src/*.h $(gtest_headers)

# Flags for gtest
test_flags := -I$(gtest_dir)/include -I$(gtest_dir)
test_lib_flags := -pthread

progs := itp2d run_tests
libs := libitp2d.a
src := $(wildcard src/*.cpp)
hdr := $(wildcard src/*.hpp)
obj := $(patsubst src/%.cpp,obj/%.o,$(src)) obj/gtest-all.o
dep := $(patsubst obj/%.o,.deps/%.o.d,$(obj)) .deps/itp2d.d .deps/run_tests.d
lib_objs := $(filter-out obj/itp2d.o obj/run_tests.o obj/test_% obj/gtest-%, $(obj))

# Make targets and rules follow

.PHONY: default check doc all lib clean depend internalchecks

default: itp2d

check: run_tests
	./run_tests

doc: README.html

lib: libitp2d.a

all: itp2d run_tests lib

clean:
	rm -f $(progs) $(libs) $(obj) $(dep)

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
	$(CXX) $(inc_flags) -MM -MT obj/$*.o $< -MF $@

.deps/gtest-all.o.d: $(gtest_srcs) scripts/depmunger.py | .deps
	$(CXX) $(inc_flags) $(test_flags) -MM -MT obj/gtest-all.o $< -MF $@

.deps/itp2d.d: src/itp2d.cpp scripts/depmunger.py | .deps
	$(CXX) $(inc_flags) -MM -MT itp2d $< | scripts/depmunger.py > $@

.deps/run_tests.d: src/run_tests.cpp scripts/depmunger.py | .deps
	$(CXX) $(inc_flags) -MM -MT run_tests $< | scripts/depmunger.py > $@

obj/test_%.o: src/test_%.cpp $(gtest_headers)| obj
	$(CXX) $(flags) $(inc_flags) $(test_flags) -c $< -o $@

obj/run_tests.o: src/run_tests.cpp $(gtest_headers) | obj
	$(CXX) $(flags) $(inc_flags) $(test_flags) -c $< -o $@

obj/gtest-all.o: $(gtest_srcs) | $(gtest_dir)
	$(CXX) $(CXXFLAGS) $(inc_flags) $(test_flags) -c $(gtest_dir)/src/gtest-all.cc -o $@

obj/%.o: src/%.cpp | obj
	$(CXX) $(flags) $(inc_flags) -c $< -o $@

itp2d: obj/itp2d.o .deps/itp2d.d | data internalchecks
	$(CXX) $(flags) $(filter %.o,$^) $(lib_flags) -o $@

run_tests: obj/run_tests.o obj/gtest-all.o .deps/run_tests.d | data internalchecks
	$(CXX) $(flags) $(filter %.o,$^) $(lib_flags) $(test_lib_flags) -o $@

libitp2d.a: $(lib_objs) | internalchecks
	$(AR) rcs $@ $^

%.html: %.md
ifdef markdown
	$(markdown) $< > $@
else
	$(info No markdown compiler specified. Cannot build HTML documentation)
endif

# include dependencies

-include $(dep)
