#!/bin/bash
################################################################
####
#### Install all variants of Petsc given loaded compiler/mpi
####
################################################################

function usage() {
    echo "Usage: $0 [ -h ] "
    echo "    [ -v (default: ${pversion} ]"
    echo "    [ -j jpar (default: ${jcount}) ]"
    echo "    [ -c : cuda build ]"
    echo "    [ -4 : include petsc/slepc4py ]"
}

cuda=0
jcount=6
p4p=0
pversion=3.22.0
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-4" ] ; then 
	p4p=1 && shift
    elif [ $1 = "-c" ] ; then 
	cuda=1 && shift
    elif [ $1 = "-j" ] ; then
	shift && jcount=$1 && shift
    elif [ $1  = "-v" ] ; then 
	shift && pversion=$1 && shift
    else 
	echo "Error: unknown option <<$1>>" && exit 1
    fi
done

somelog=some.log
rm -f ${somelog} && touch ${somelog}

( \
echo "================================================================" && \
echo "" && \
echo "    Installation petsc variants for:" && \
echo "    compiler=${TACC_FAMILY_COMPILER}/${TACC_FAMILY_COMPILER_VERSION} mpi=${TACC_FAMILY_MPI}/${TACC_FAMILY_MPI_VERSION}" && \
echo "" && \
echo "================================================================" \
) | tee -a ${somelog}      

set -e
archs=archs-${pversion}
rm -f $archs && touch $archs

#
# install the arch list that we test in all_local/tacc_tests.sh
#
for arch in  $( cat test_versions.txt ) ; do
    case "${arch}" in ( *debug* ) export DEBUG=1 ;; ( * ) export DEBUG=0 ;; esac
    case "${arch}" in ( *i64 ) export INT=64 ;; ( * ) export INT=32 ;; esac
    case "${arch}" in ( *single* ) export PRECISION=single ;; ( * ) export PRECISION=double ;; esac
    case "${arch}" in ( *complex* ) export SCALAR=complex ;; ( * ) export SCALAR=real ;; esac
    arch=$( make --no-print-directory petscshortarch )
    ( echo && echo "Installing big for arch=${arch}" && echo ) | tee -a ${somelog}
    if [ -z "${arch}" ] ; then 
	echo vanilla >> ${archs}
    else
	echo ${arch} >> ${archs}
    fi
    ./install_big.sh -j ${jcount} -v ${pversion} \
	$( if [ "${arch}" = "f08" ] ; then echo "-8" ; fi ) \
	$( if [ "${arch}" = "nohdf5" ] ; then echo "-5" ; fi ) \
	$( if [ ${p4p} -eq 1 ] ; then echo "-4" ; fi ) \
	$( if [ ${cuda} -eq 1 ] ; then echo "-c" ; fi ) \
	| tee -a ${somelog}
done


( echo && echo "done archs: $( cat ${archs} | tr '\n' ' ' )" && echo ) | tee -a ${somelog}
