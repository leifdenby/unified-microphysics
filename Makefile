objects := $(patsubst %.f90,%.o, $(wildcard mphys_*.f90))

all: $(objects) libunified_microphysics.a
	echo $(objects)

libunified_microphysics.a: $(objects)
	python generate_wrapper.py
	gfortran -c microphysics_register.f90
	ar rc libunified_microphysics.a microphysics_register.o $(objects) *.mod

microphysics_common.o: microphysics_common.f90
	gfortran -c microphysics_common.f90

test: all libunified_microphysics.a
	#gfortran -l unified_microphysics -L. test.f90 -o test
	gfortran -g microphysics_register.o -l unified_microphysics -L. test.f90 -o test
	#gfortran mphys_dummy.o microphysics_register.o test.f90 -o test
	./test

$(objects): $(wildcard mphys_*.f90) microphysics_common.o
	gfortran -c $(patsubst %.o, %.f90, $@)

clean:
	rm *.o *.mod microphysics_register.f90 test lib*.a
