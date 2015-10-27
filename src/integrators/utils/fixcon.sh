#!/bin/bash
#
gfortran -c fixcon.f90
if [ $? -ne 0 ]; then
  echo "Errors while compiling fixcon.f90"
  exit
fi
#
gfortran fixcon.o
if [ $? -ne 0 ]; then
  echo "Errors while loading fixcon.o"
  exit
fi
rm fixcon.o
#
mv a.out ~/bin/$ARCH/fixcon
#
echo "Executable installed as ~/bin/$ARCH/fixcon."
