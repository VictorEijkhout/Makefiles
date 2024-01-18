#!/bin/bash

jcount=6
pversion=3.20.0
function usage() {
    echo "Usage: $0 [ -h ] [ -j 123 ] [ -v (default: ${pversion} ]"
    echo "    [ -4 : skip petsc/slepc4py ]"
    echo "    [ -f : skip fortran ]"
}

p4p=1
fort=1
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-4" ] ; then 
	p4p=0 && shift
    elif [ $1 = "-f" ] ; then 
	fort=0 && shift
    elif [ $1 = "-j" ] ; then
	shift && jcount=$1 && shift
    elif [ $1  = "-v" ] ; then 
	shift && pversion=$1 && shift
    else 
	echo "Error: unknown option <<$1>>" && exit 1
    fi
done

if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then
    module load mkl
fi
module load eigen
module load fftw3
#module load libceed
module load phdf5/1.14

make --no-print-directory allinstall \
     PACKAGEVERSION=${pversion} JCOUNT=${jcount} \
     FORTRAN=${fort} PETSC4PY=${p4p} SLEPC4PY=${p4p}

