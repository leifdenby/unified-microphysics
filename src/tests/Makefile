# Notes
#
# - `-fPIC` flag is necessary to create library that can be externally linked, from e.g. python
# - library has to go at end of compiler command when linking against static library
#     i.e. not
#     	$(GC) $(FFLAGS) microphysics_tests.F90 -o microphysics_tests -lunified_microphysics -L.
#     instead
#     	$(GC) $(FFLAGS) -lunified_microphysics -L. microphysics_tests.F90 -o microphysics_tests

FFLAGS :=-g -O0 -fPIC
GC:=gfortran  # must use gfortran until I find out how to make f2py use ifort

all: python
	#$(GC) -l unified_microphysics -L. test.F90 -o test
	$(GC) $(FFLAGS) microphysics_tests.F90 -o microphysics_tests -lunified_microphysics -L../base -I../base
	#gfortran mphys_dummy.o microphysics_register.o test.F90 -o test
	./microphysics_tests
	env PYTHONPATH=../../pylib/:$(PYTHONPATH) py.test ../../pylib/unified_microphysics/tests -v

python:
	$(MAKE) -C ../wrappers python

clean:
	rm microphysics_tests
	rm *.o *.mod lib*.a
