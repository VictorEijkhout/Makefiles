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
module load phdf5/1.14

module load cuda/12

export biglog=big_install.log
for d in 0 1 ; do 
    make --no-print-directory biginstall EXT= JCOUNT=${jcount} \
	 DEBUG=$d \
	 EIGEN=1 FFTW3=1 HDF5=1 HYPRE=0 \
	 CUDA=0 FORTRAN=0 \
	 PETSC4PY=1 SLEPC4PY=1 \
	 2>&1 | tee -a ${biglog}
    make --no-print-directory biginstall EXT= JCOUNT=${jcount} \
	 DEBUG=$d \
	 CUDA=1 FORTRAN=0 \
	 2>&1 | tee -a ${biglog}
done

echo && echo "See: ${biglog}" && echo
    
