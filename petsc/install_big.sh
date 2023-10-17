#!/bin/bash

jcount=6
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-j" ] ; then
	shift && jcount=$1 && shift
    else 
	echo "Error: unknown option <<$1>>" && exit 1
    fi
done

if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then
    module load mkl
fi
module load eigen
module load fftw3
module load hypre
## module load libceed
module load phdf5/1.14
echo "hdf5 lib: $TACC_PHDF5_LIB"

make --no-print-directory biginstall EXT= JCOUNT=${jcount} \
    EIGEN=1 FFTW3=1 HDF5=1 HYPRE=0 \
    CUDA=0 FORTRAN=0 \
    PETSC4PY=1 SLEPC4PY=1
