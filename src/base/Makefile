model_sources = $(wildcard ../models/mphys_*.F90)
objects := $(patsubst %.F90,%.o, $(model_sources))

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

all: $(objects) libunified_microphysics.a

microphysics_constants.o: microphysics_constants.F90
	$(GC) $(FFLAGS) -c microphysics_constants.F90

microphysics_register.o: microphysics_constants.o microphysics_register.F90
	$(GC) $(FFLAGS) -c microphysics_register.F90

microphysics_common.o: microphysics_register.o microphysics_common.F90
	$(GC) $(FFLAGS) -c microphysics_common.F90

microphysics_initialisation.o: $(model_sources)
	python generate_wrapper.py
	$(GC) $(FFLAGS) -c microphysics_initialisation.F90

# module for each microphysics implementation
$(objects): $(wildcard mphys_*.F90) microphysics_common.o microphysics_register.o
	$(GC) $(FFLAGS) -c $(patsubst %.o, %.F90, $@)

libunified_microphysics.a: $(objects) microphysics_register.o microphysics_initialisation.o
	ar crs libunified_microphysics.a microphysics_register.o microphysics_common.o microphysics_initialisation.o $(subst ../models/,,$(objects))

test: all libunified_microphysics.a
	#$(GC) -l unified_microphysics -L. test.F90 -o test
	$(GC) $(FFLAGS) microphysics_tests.F90 -o microphysics_tests -lunified_microphysics -L.
	#gfortran mphys_dummy.o microphysics_register.o test.F90 -o test
	./microphysics_tests

clean:
	rm *.o *.mod lib*.a
	rm microphysics_initialisation.F90
