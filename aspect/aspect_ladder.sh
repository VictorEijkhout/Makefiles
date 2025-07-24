#!/bin/bash

setx=0
jcount=4
list=
packages=0

phdf5_dir=hdf5
phdf5_tgt=par
PACKAGEOPTIONS_hdf5="HDFFORTRAN=OFF"

petsc_commandline="./install_big.sh -f"
petsc_version=3.23
petsc_full_version=3.23.2

# these are module names
ladder="\
    cmake \
    petsc,${petsc_version} \
    zlib \
    p4est \
    boost \
    pcre2 \
    bison \
    swig \
    metis \
    phdf5 \
    autoconf,2.72 \
    automake,1.16.5 \
    pnetcdf \
    trilinos \
    sundials \
    dealii \
    aspect \
    "
#     gklib is maybe no needed

source ../ladder.sh

