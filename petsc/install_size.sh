#!/bin/bash

jcount=6
pversion=
function usage() {
    echo "Usage: $0 [ -h ] [ -v : leave null for default ]"
    echo "    [ -j jpar (default: ${jcount}) ]"
    echo "    [ -e customext ]"
    echo "    [ -4 : include petsc/slepc4py ] [ -5 : skip hdf5 ]"
    echo "    [ -f : skip fortran ] [ -8 : f08 support ]"
    echo "    [ -k : kokkos support ] [ -u : cuda build ]"
    echo "    environment: DEBUG INT PRECISION SCALAR"
}

echo && echo "Starting ${size} installation"
echo " .. compiler=${TACC_FAMILY_COMPILER}/${TACC_FAMILY_COMPILER_VERSION}"

fortran=90
hdf5=1
p4p=0
cuda=0
kokkos=0
customext=
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-4" ] ; then 
	echo " .. enabling petsc4py/slepc4py"
	p4p=1 && shift
    elif [ $1 = "-5" ] ; then 
	echo " .. disabling hdf5"
	hdf5=0 && shift
    elif [ $1 = "-8" ] ; then 
	echo " .. enable f08"
	fortran=08 && shift
    elif [ $1 = "-e" ] ; then 
	shift && customext=$1 && shift
	echo " .. custom extension: $customext"
    elif [ $1 = "-f" ] ; then 
	echo " .. disable fortran"
	fortran=0 && shift
    elif [ $1 = "-j" ] ; then
	shift && jcount=$1 && shift
    elif [ $1 = "-k" ] ; then
	kokkos=1 && shift
    elif [ $1 = "-u" ] ; then 
	echo " .. using CUDA"
	cuda=1 && shift
    elif [ $1  = "-v" ] ; then 
	shift && pversion=$1 && shift
    else 
	echo "Error: unknown option <<$1>>" && exit 1
    fi
done

if [ $cuda -eq 1 ] ; then
    echo " .. disabling Fortran for Cuda"
    fortran=0
fi

echo " .. make thread count $jcount"

set -e
case "${TACC_SYSTEM}" in
  ( "vista" )
    module load nvpl ;;
  ( * )
    if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then
      module load mkl
    fi ;;
esac

if [ "${cuda}" = "1" ] ; then 
    gcc_version=${TACC_FAMILY_COMPILER_VERSION}
    gcc_version=${gcc_version%%.*}
    echo "Loading CUDA module"
    module load cuda/12.8
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

export sizelog=install_${size}$( if [ ! -z "${customext}" ] ; then echo "-${customext}" ; fi ).log
rm -f $sizelog
cmdline="\
make --no-print-directory biginstall JCOUNT=${jcount} \
    $( if [ ! -z "${pversion}" ] ; then echo PACKAGEVERSION=${pversion} ; fi ) \
    EXT=${EXTENSION} \
    $( if [ ! -z "${customext}" ] ; then echo CUSTOMEXT=${customext} ; fi ) \
    PKGSTART=1 $( cat ${size}.sh ) PKGSTOP=1 \
    CUDA=${cuda} FORTRAN=${fortran} \
    $( if [ "${kokkos}" = "1" -o "${KOKKOS}" = "1" ] ; then echo KOKKOS=1 ; fi ) \
    PETSC4PY=${p4p} SLEPC4PY=${p4p} \
"
( echo "At $(date)" && echo "cmdline: $cmdline" ) | tee -a ${sizelog}
set -e
set -o pipefail 
eval $cmdline 2>&1 | tee -a ${sizelog}
# if [ $retcode -gt 0 ] ; then
#     echo "Make failed for EXT=${EXTENSION}" 
#     echo && echo "See: ${sizelog}" && echo
#     exit 1
# else
    echo && echo "See: ${sizelog}" && echo
##fi

