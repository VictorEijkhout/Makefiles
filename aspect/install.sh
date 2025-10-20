#!/bin/bash

set -e

petsc=1
while [ $# -gt 0 ] ; do
    if [ $1 = "-p" ] ; then
	petsc=0 && shift
    fi
done


if [ $petsc = "1" ] ; then 
    cd ../petsc
    ./install_big.sh -f
fi
module load petsc/3.24

module load boost
module load swig
module load phdf5
module load pnetcdf

cd ../trilinos
make default_install
module load trilinos

cd ../dealii
module load p4est sundials
make default_install
