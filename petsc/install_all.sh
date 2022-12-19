#!/bin/bash

version=3.18.2
if [ $# -eq 1 ] ; then 
  version=$1
fi

if [ "${TACC_SYSTEM}" != "macbookair" ] ; then
    options="P4P=1 KOKKOS=1 FORTRAN=0"
    module load cuda
    for d in 1 0 ; do
	make install DEBUG=$d PACKAGEVERSION=$version \
	     $options CUDA=1
    done
    module unload cuda
fi

if [ "${TACC_SYSTEM}" = "ls6" -o "${TACC_SYSTEM}" = "macbookair" ] ; then
    fortran=0
else
    fortran=1
fi
options="HDF5=1 METIS=0 P4P=1 KOKKOS=0 SLEPC=1 FORTRAN=$fortran"
for d in 1 0 ; do
    # make with default official version
    make install DEBUG=$d PACKAGEVERSION=$version \
	$options
    # make latest git pull
    # make configure install DEBUG=$d $options PACKAGEVERSION=git
done

