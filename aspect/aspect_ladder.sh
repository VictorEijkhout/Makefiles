#!/bin/bash

setx=0
jcount=4
list=
packages=0

phdf5_dir=hdf5
phdf5_tgt=par
PACKAGEOPTIONS_hdf5="HDFFORTRAN=OFF"

petsc_version=3.23
petsc_full_version=3.23.1

trilinosversion=14.4.0
# these are module names
ladder="\
    zlib \
    petsc,3.23 \
    p4est \
    boost \
    pcre2 \
    bison \
    swig \
    gklib \
    metis \
    phdf5 \
    pnetcdf \
    trilinos \
    dealii \
    aspect \
    "

source ../ladder.sh

