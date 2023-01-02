#!/bin/bash

version=3.18.2
if [ $# -eq 1 ] ; then 
  version=$1
fi

if [ "${TACC_SYSTEM}" = "ls6" -o "${TACC_SYSTEM}" = "macbookair" ] ; then
    fortran=0
else
    fortran=1
fi
if [ "${TACC_SYSTEM}" = "stampede2" ] ; then 
    p4p=0
else
    p4p=1
fi
slepc=0

##
## regular versions
##
options="HDF5=1 METIS=0 P4P=${p4p} KOKKOS=0 SLEPC=${slepc} FORTRAN=${fortran}"
for d in 0 1 ; do
    # make with default official version
    make install DEBUG=$d PACKAGEVERSION=$version \
	$options
    exit
done

##
## cuda versions
##
if [ "${TACC_SYSTEM}" = "frontera" -o "${TACC_SYSTEM}" = "ls6" ] ; then
    options="P4P=${p4p} KOKKOS=1 FORTRAN=0"
    module load cuda
    for d in 1 0 ; do
	make install DEBUG=$d PACKAGEVERSION=$version \
	     $options CUDA=1
    done
    module unload cuda
fi
