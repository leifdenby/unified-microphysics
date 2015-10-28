integrator_sources = $(wildcard _*.f90)
objects := $(patsubst %.f90,%.o, $(integrator_sources))

F90?=gfortran
FFLAGS :=-g -O0 -fPIC

all: $(objects)
	ar crs libmicrophysics_integrators.a $(objects)

$(objects): $(integrator_sources)
	$(F90) $(FFLAGS) -c $<

test: test.f90 $(objects)
	$(F90) _rkf45.o test.f90 -o test
	./test

clean:
	rm -f *.o *.mod *.a
