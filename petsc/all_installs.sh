#!/bin/bash

version=3.18.2
if [ $# -eq 1 ] ; then 
  version=$1
fi

if [ "${TACC_SYSTEM}" != "macbookair" ] ; then
    options="P4P=1 KOKKOS=1"
    module load cuda
    for d in 1 0 ; do
	make install DEBUG=$d PACKAGEVERSION=$version \
	     $options CUDA=1
    done
    module unload cuda
fi

options="HDF5=1 METIS=1 P4P=1 KOKKOS=1 SLEPC=1 FFTW=1"
for d in 1 0 ; do
    # make with default official version
    make install DEBUG=$d PACKAGEVERSION=$version \
	$options
    # make latest git pull
    # make configure install DEBUG=$d $options PACKAGEVERSION=git
done

