#!/bin/bash

setx=0
jcount=4
list=
packages=0

ladder="\
    zlib,1.2.13 \
    phdf5,1.14.3,hdf5,par \
    parallelnetcdf,git \
    pnetcdf,4.9.2,netcdf,par \
    boost,1.83.0 \
    trilinos,15.0.0 \
    peridigm,git \
    "
