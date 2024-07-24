#!/bin/bash

setx=0
jcount=4
list=
packages=0

hdf5_target=par
PACKAGEOPTIONS_hdf5="HDFFORTRAN=OFF"

petsc_version=3.20
petsc_full_version=3.20.5

trilinosversion=14.4.0
ladder="\
    zlib,1.2.13 \
    petsc,${petsc_version} \
    p4est,2.8.6 \
    boost,1.83.0 \
    pcre2,git \
    swig,4.1.1 \
    hdf5,1.14 \
    netcdf,4.9.2 \
    gklib,git \
    metis,5.1.0.3 \
    trilinos,${trilinosversion} \
    dealii,9.5.1 \
    aspect,2.5.0 \
    "

source ../ladder.sh

