#!/bin/bash

setx=0
jcount=4
list=
packages=0

hdf5_target=par
PACKAGEOPTIONS_hdf5="HDFFORTRAN=OFF"

petsc_version=3.23
petsc_full_version=3.23.1

trilinosversion=14.4.0
ladder="\
    zlib \
    petsc,3.23 \
    p4est \
    boost \
    pcre2 \
    bison \
    swig \
    hdf5 \
    netcdf \
    gklib \
    metis \
    trilinos \
    dealii \
    aspect \
    "

source ../ladder.sh

