#!/bin/bash

jcount=6
pversion=3.20.4
function usage() {
    echo "Usage: $0 [ -h ] [ -v (default: ${pversion} ]"
    echo "    [ -j jpar (default: ${jcount}) ]"
    echo "    [ -c : cuda build ]"
    echo "    [ -e customext ]"
    echo "    [ -4 : skip petsc/slepc4py ] [ -5 : skip hdf5 ] [ -8 : f08 suupport ]"
    echo "    environment: DEBUG INT PRECISION SCALAR"
}

echo && echo "Starting big installation"
echo " .. compiler=${TACC_FAMILY_COMPILER}/${TACC_FAMILY_COMPILER_VERSION}"

fortran=90
hdf5=1
p4p=1
cuda=0
customext=
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-4" ] ; then 
	echo " .. disabling petsc4py/slepc4py"
	p4p=0 && shift
    elif [ $1 = "-5" ] ; then 
	echo " .. disabling hdf5"
	hdf5=0 && shift
    elif [ $1 = "-8" ] ; then 
	echo " .. enable f08"
	fortran=08 && shift
    elif [ $1 = "-c" ] ; then 
	echo " .. using CUDA"
	cuda=1 && shift
    elif [ $1 = "-e" ] ; then 
	shift && customext=$1 && shift
	echo " .. custom extension: $customext"
    elif [ $1 = "-j" ] ; then
	shift && jcount=$1 && shift
    elif [ $1  = "-v" ] ; then 
	shift && pversion=$1 && shift
    else 
	echo "Error: unknown option <<$1>>" && exit 1
    fi
done
echo " .. make thread count $jcount"

set -e
if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then
    module load mkl
fi

module load eigen
module load fftw3
if [ "${INT}" = "64" ] ; then 
    module load hypre/2.30.0-i64
else
    module load hypre/2.30.0
fi
## module load phdf5
# /1.14

if [ "${cuda}" = "1" ] ; then 
    cversion=${TACC_FAMILY_COMPILER_VERSION}
    cversion=${cversion%%.*}
    if [ ${cversion} -gt 12 ] ; then 
	echo "ERROR can not deal with gcc > 12" && exit 1
    fi
    module load cuda/12
fi

EXTENSION=
if [ "${hdf5}" = "0" ] ; then
    EXTENSION=${EXTENSION}nohdf5
fi
if [ "${SCALAR}" = "complex" ] ; then 
    EXTENSION=${EXTENSION}complex
fi
if [ "${PRECISION}" = "single" ] ; then 
    EXTENSION=${EXTENSION}single
fi
if [ "${INT}" = "64" ] ; then 
    EXTENSION=${EXTENSION}i64
fi
if [ "${fortran}" = "08" ] ; then 
    EXTENSION=${EXTENSION}f08
fi
if [ "${DEBUG}" = "1" ] ; then
    EXTENSION=${EXTENSION}debug
fi
if [ "${TACC_FAMILY_COMPILER}" = "intel" ] ; then
    echo " .. Note: Intel disabling CHACO"
    CHACO=0
else
    CHACO=1
fi

export biglog=install_big$( if [ ! -z "${customext}" ] ; then echo "-${customext}" ; fi ).log
rm -f $biglog
cmdline="\
make --no-print-directory biginstall JCOUNT=${jcount} PACKAGEVERSION=${pversion} \
    EXT=${EXTENSION} \
    $( if [ ! -z "${customext}" ] ; then echo CUSTOMEXT=${customext} ; fi ) \
    AMGX=1 CHACO=${CHACO} EIGEN=1 FFTW3=1 HDF5=${hdf5} HYPRE=1 MUMPS=1 METIS=1 \
    PARMETIS=1 PTSCOTCH=1 \
    CUDA=${cuda} FORTRAN=${fortran} \
    PETSC4PY=${p4p} SLEPC4PY=${p4p} \
"
echo "At $(date) cmdline: $cmdline" | tee -a ${biglog}
set -e
set -o pipefail 
eval $cmdline 2>&1 | tee -a ${biglog}
# if [ $retcode -gt 0 ] ; then
#     echo "Make failed for EXT=${EXTENSION}" 
#     echo && echo "See: ${biglog}" && echo
#     exit 1
# else
    echo && echo "See: ${biglog}" && echo
##fi

