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

#FFLAGS :=-g -O0 -fPIC
FFLAGS:=-O3 -fPIC
FC?=gfortran  # must use gfortran until I find out how to make f2py use ifort

all: $(objects) libunified_microphysics.a

ifeq ($(FC),mpif90)
FFLAGS += -DMPI
endif

../integrators/libmicrophysics_integrators.a:
	FC=$(FC) $(MAKE) -C ../integrators

microphysics_constants.o: microphysics_constants.F90
	$(FC) $(FFLAGS) -c microphysics_constants.F90

microphysics_register.o: microphysics_constants.o
	$(FC) $(FFLAGS) -c microphysics_register.F90

microphysics_common.o: microphysics_register.o microphysics_common.F90
	$(FC) $(FFLAGS) -c microphysics_common.F90

integrators.o: integrators.F90
	$(FC) $(FFLAGS) -c integrators.F90

microphysics_integration.o: microphysics_register.o microphysics_common.F90 ../integrators/libmicrophysics_integrators.a integrators.o
	$(FC) $(FFLAGS) -c microphysics_integration.F90 -L../integrators -lmicrophysics_integrators -I../integrators

microphysics_initialisation.o: $(model_sources) microphysics_initialisation.tpl.F90 microphysics_integration.o
	python generate_wrapper.py
	$(FC) $(FFLAGS) -c microphysics_initialisation.F90

# module for each microphysics implementation
$(objects): $(wildcard mphys_*.F90) microphysics_common.o microphysics_register.o microphysics_constants.o
	$(FC) $(FFLAGS) -c $(patsubst %.o, %.F90, $@)

libunified_microphysics.a: $(objects) microphysics_register.o microphysics_initialisation.o microphysics_integration.o
	ar crs libunified_microphysics.a microphysics_register.o microphysics_common.o microphysics_initialisation.o microphysics_integration.o integrators.o $(subst ../models/,,$(objects))
	mkdir -p ../../lib
	cp libunified_microphysics.a ../../lib
	mkdir -p ../../include
	cp *.mod ../../include

clean:
	rm -f *.o *.mod lib*.a
	rm -f microphysics_initialisation.F90
	rm -rf ../../include
	rm -rf ../../lib
