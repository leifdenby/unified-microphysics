_Aim_: microphysics module written in Fortran which can be compiled stand-alone or integrated into ATHAM, CCFM, KiD, etc.

# Preliminary spec

Internally the state-vector `y` represents all the scalar fields relevant to a
specific microphysics model in a single computational volume. Given the initial
state `y0` at a time `t0` it is the aim of the microphysics framework to
calculate the new state at time `t1` later. This integration can either be
performed isobarically or isometrically.

The microphysics framework consists of the following components:

1. A number of microphysics "models" which each serve to calculate the time
   derivative of the state variable (dQ/dt). These are in [src/models/](src/models/)
2. A central "register" which handles the booking keeping with respect to
   picking at runtime a specific microphyiscs model and well as an abstraction
   which guarantees physical consistency when integrating the selected model.
   Resides in [src/base/](src/base/)
3. A collection of numerical integration routines, implemented independently of
   the underlying physics.  Currently only the Runge-Kutta 3-4 is implemented.
   These are in [src/base/integrators.f90](src/base/integrators.f90) but will move to [src/integrators/](src/integrators/)
   eventually.
4. A number of utility functions for calculating thermodynamic variables and
   common parameterisations. These are in [src/base/microphysics_common.f90](src/base/microphysics_common.f90).

# Using the framework

There are currently two main ways of using the microphysics framework, either
as an external library within a Fortran90 program or as a module from Python.

## In an external Fortran90 program

1. Compile (set compiler with `FC` environment variable)

        FC=ifort make

2. Write wrapper to initiate the microphysics module and to call the
   integration in every timestep of the host model. See the ATHAM wrapper in
   [https://github.com/leifdenby/unified-microphysics](src/wrappers/atham_microphysics.f90) for inspiration.

3. When compiling use header files in `include/` and link against library in `lib/`, i.e.

        export UM_PATH=/home/lcd33/git-repos/unified-microphysics
        ifort external_program.f90 -I$UM_PATH/include -L$UM_PATH/lib -lunified_microphysics

## From Python

1. Compile Python wrapper
    
        FC=ifort make python

2. Add `pylib/` to `$PYTHONPATH` environment variable, e.g.

        export PYTHONPATH=$PYTHONPATH:/home/lcd33/git-repos/unified-microphysics/pylib

3. Import `unified_microphysics` module in python and use. All tests are
   currently written in Python so these (in
   [pylib/unified_microphysics/tests/](pylib/unified_microphysics/tests/)) may serve as inspiration.


## Running tests

All tests are currently written in Python and require `py.test` to run so this
will first need installing (eg with `pip install pytest`) and then tests can be
run with

        FC=ifort make tests


# Resources

On Makefiles:

- http://stackoverflow.com/questions/8937500/how-to-generate-list-of-make-targets-automatically-by-globbing-subdirectories
- http://www.gnu.org/software/make/manual/make.html#Wildcard-Function
- http://stackoverflow.com/questions/1633527/wildcard-targets-in-a-makefile

On Fortran module compilation:

- http://docs.oracle.com/cd/E19205-01/819-5263/aevog/index.html
