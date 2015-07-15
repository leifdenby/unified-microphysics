objects := $(patsubst %.f90,%.o, $(wildcard mphys_*.f90))

all: $(objects)
	echo $(objects)
	python generate_wrapper.py
	gfortran $(objects) microphysics_register.f90 test.f90 -o test

test: all
	./test

$(objects): $(wildcard mphys_*.f90)
	gfortran -c $(patsubst %.o, %.f90, $@)

clean:
	rm *.o *.mod microphysics_register.f90
