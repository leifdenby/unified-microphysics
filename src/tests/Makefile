# Notes
#
# - `-fPIC` flag is necessary to create library that can be externally linked, from e.g. python
# - library has to go at end of compiler command when linking against static library
#     i.e. not
#     	$(GC) $(FFLAGS) microphysics_tests.F90 -o microphysics_tests -lunified_microphysics -L.
#     instead
#     	$(GC) $(FFLAGS) -lunified_microphysics -L. microphysics_tests.F90 -o microphysics_tests

FFLAGS :=-g -O0 -fPIC
FC:=ifort

all: python
	$(FC) $(FFLAGS) microphysics_tests.F90 -o microphysics_tests -lunified_microphysics -L../base -I../base
	./microphysics_tests
	env PYTHONPATH=../../pylib/:$(PYTHONPATH) py.test ../../pylib/unified_microphysics/tests -v

python:
	$(MAKE) -C ../wrappers python

clean:
	rm microphysics_tests
	rm *.o *.mod lib*.a
