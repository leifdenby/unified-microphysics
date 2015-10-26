objects := $(patsubst %.F90,%.o, $(wildcard mphys_*.F90))

# Notes
#
# - `-fPIC` flag is necessary to create library that can be externally linked, from e.g. python
# - library has to go at end of compiler command when linking against static library
#     i.e. not
#     	$(GC) $(FFLAGS) microphysics_tests.F90 -o microphysics_tests -lunified_microphysics -L.
#     instead
#     	$(GC) $(FFLAGS) -lunified_microphysics -L. microphysics_tests.F90 -o microphysics_tests

FFLAGS :=-g -O0 -fPIC
FC?=gfortran  # must use gfortran until I find out how to make f2py use ifort

all: base

base:
	FC=$(FC) $(MAKE) -C src/base

python: base
	FC=$(FC) $(MAKE) -C src/wrappers python

tests: base
	FC=$(FC) $(MAKE) -C src/tests

clean:
	FC=$(FC) $(MAKE) -C src/base clean
	FC=$(FC) $(MAKE) -C src/wrappers clean
	rm -f lib/*
	rm -f include/*

