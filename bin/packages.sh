#!/bin/bash
################################################################
####
#### List my favorite packages
####
################################################################

cpplibs="catch2 cxxopts eigen fmtlib mdspan msgsl rangev3 sfml"
install="cmake meson ninja scons"
scicomp="adios2 arpack aspect fftw3 parmetis petsc siesta suitesparse trilinos zoltan"
storage="hdf5 phdf5 netcdf parallelnetcdf pnetcdf silo zlib"
utils="autoconf automake bison boost pcre2 swig"

for p in \
    ${cpplibs} ${install} ${scicomp} ${storage} ${utils} 
    do
    echo $p
done
	 
