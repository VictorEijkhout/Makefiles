export CC=icx
export CFLAGS="-O2 -ansi -std=c99 -fopenmp -fPIC"
export CCOPTIONS="-O2 -ansi -std=c99 -fopenmp -fPIC"

export CXX=icpx
export CXXFLAGS="-O2 -ansi -std=c99 -fopenmp -fPIC"
export CPPFLAGS="-DNDEBUG"

export FC=ifx
export FFLAGS='-O2 -fPIC -fopenmp -fallow-argument-mismatch -fallow-invalid-boz'
export FCOPTIONS='-O2 -fPIC -fopenmp -fallow-argument-mismatch -fallow-invalid-boz'

export F90=ifx
export F90FLAGS="-O2 -fPIC -fopenmp -fallow-argument-mismatch"

