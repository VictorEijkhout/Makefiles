#!/bin/bash

jcount=6
pversion=3.20.2
function usage() {
    echo "Usage: $0 [ -h ] [ -v (default: ${pversion} ]"
    echo "    [ -j jpar (default: ${jcount}) ]"
    echo "    [ -c : cuda build ]"
    echo "    [ -e customext ]"
    echo "    [ -4 : skip petsc/slepc4py ]"
}

p4p=1
cuda=0
customext=
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-4" ] ; then 
	p4p=0 && shift
    elif [ $1 = "-c" ] ; then 
	cuda=1 && shift
    elif [ $1 = "-e" ] ; then 
	shift && customext=$1 && shift
    elif [ $1 = "-j" ] ; then
	shift && jcount=$1 && shift
    elif [ $1  = "-v" ] ; then 
	shift && pversion=$1 && shift
    else 
	echo "Error: unknown option <<$1>>" && exit 1
    fi
done

set -e
if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then
    module load mkl
fi

module load eigen
module load fftw3
if [ "${INT}" = "64" ] ; then 
    module load hypre/2.29.0-i64
else
    module load hypre/2.29.0
fi
module load phdf5/1.14

if [ "${cuda}" = "1" ] ; then 
    cversion=${TACC_FAMILY_COMPILER_VERSION}
    cversion=${cversion%%.*}
    if [ ${cversion} -gt 12 ] ; then 
	echo "ERROR can not deal with gcc > 12" && exit 1
    fi
    module load cuda/12
fi

export biglog=install_big$( if [ ! -z "${customext}" ] ; then echo "-${customext}" ; fi ).log
rm -f $biglog
EXTENSION=
if [ "${SCALAR}" = "complex" ] ; then 
    EXTENSION=${EXTENSION}complex
fi
if [ "${PRECISION}" = "single" ] ; then 
    EXTENSION=${EXTENSION}single
fi
if [ "${INT}" = "64" ] ; then 
    EXTENSION=${EXTENSION}i64
fi
if [ "${DEBUG}" = "1" ] ; then
    EXTENSION=${EXTENSION}debug
fi
cmdline="\
make --no-print-directory biginstall JCOUNT=${jcount} PACKAGEVERSION=${pversion} \
    EXT=${EXTENSION} \
    $( if [ ! -z "${customext}" ] ; then echo CUSTOMEXT=${customext} ; fi ) \
    AMGX=1 CHACO=1 EIGEN=1 FFTW3=1 HDF5=1 HYPRE=1 PARMETIS=1 \
    CUDA=${cuda} FORTRAN=1 \
    PETSC4PY=${p4p} SLEPC4PY=${p4p} \
"
echo "cmdline: $cmdline" | tee -a ${biglog}
eval $cmdline  2>&1 | tee -a ${biglog}

echo && echo "See: ${biglog}" && echo
    
