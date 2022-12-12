#!/bin/bash

################################################################
####
#### installation of complex petsc versions for deal
####
################################################################

version=3.18.1
if [ $# -eq 1 ] ; then 
  version=$1
fi

if [ "${TACC_SYSTEM}" = "ls6" ] ; then
    fortran=0
else
    fortran=1
fi
options="SCALAR=complex SLEPC=1 FORTRAN=$fortran"
module unload cuda
for d in 1 0 ; do
    make install DEBUG=$d PACKAGEVERSION=$version \
	$options
done

