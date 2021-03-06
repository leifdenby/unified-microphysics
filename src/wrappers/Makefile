models := $(wildcard ../models/mphys_*.F90)
objects := $(patsubst %.F90,%.o, $(models))

base_src = ../base/microphysics_constants.F90 ../base/microphysics_common.F90 ../base/microphysics_register.F90

# Notes
#
# - `-fPIC` flag is necessary to create library that can be externally linked, from e.g. python
# - library has to go at end of compiler command when linking against static library
#     i.e. not
#     	$(GC) $(FFLAGS) microphysics_tests.F90 -o microphysics_tests -lunified_microphysics -L.
#     instead
#     	$(GC) $(FFLAGS) -lunified_microphysics -L. microphysics_tests.F90 -o microphysics_tests

FFLAGS :=-g -O0 -fPIC
FC?=ifort

IFLAGS:=-I../base -I../integrators
ifeq ($(FC),mpif90)
	FCOMPILER=intelem
else ifeq ($(FC),ifort)
	FCOMPILER=intelem
else ifeq ($(FC),gfortran)
	FCOMPILER=gfortran
endif

unified_microphysics_fortran.so: microphysics_pylib.F90 $(models) $(base_src)
	#f2py -c -m unified_microphysics microphysics_constants.F90 microphysics_common.F90 microphysics_register.F90  #microphysics_initialisation.o $(objects)
	#f2py -L../base -lunified_microphysics -I../base --no-lower -c -m unified_microphysics microphysics_pylib.F90 skip: cv_mixture cp_mixture array_append_real array_append_int :

	# TODO: Can't import `calc_dqdt` right now because the compilation requires the number of hydrometeors set, maybe we can define a constant somewhere?
	#f2py -L../base -lunified_microphysics -I../base --no-lower -c -m unified_microphysics_fortran $(base_src) $(models) microphysics_pylib.F90 skip: cv_mixture cp_mixture array_append_real array_append_int calc_dqdt f advance :
	f2py --fcompiler=$(FCOMPILER) -L../base -lunified_microphysics $(IFLAGS) --no-lower -c -m unified_microphysics_fortran $(base_src) $(models) microphysics_pylib.F90 skip: integrate cv_mixture cp_mixture cp_gas cv_gas array_append_real array_append_int calc_dqdt f advance rkf34_original dydt :
	cp unified_microphysics_fortran.so ../../pylib/unified_microphysics/

python: unified_microphysics_fortran.so

clean:
	rm -f *.o *.mod lib*.a *.so
