#!/bin/bash

jcount=6
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-j" ] ; then
	shift && jcount=$1 && shift
    fi
done

module load mkl
module load eigen
module load fftw3
## module load libceed
module load phdf5/1.14

make --no-print-directory biginstall JCOUNT=${jcount} \
     PETSC4PY=1 SLEPC4PY=1

