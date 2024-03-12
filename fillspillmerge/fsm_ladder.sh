#!/bin/bash

setx=0
jcount=4
list=
packages=0

ladder="\
    zlib,1.2.13 \
    pcre2,10.42 \
    swig,4.1.1 \
    jsonc,git \
    sqlite,3.43.0 \
    proj,9.3.1 \
    gdal,3.7.0 \
    fillspillmerge,git \
    "

source ../ladder.sh
