################################################################
####
#### Makefile for Fftw3 installation
####
#### double variant
####
################################################################

PACKAGE = fftw3
ABOUT = Fastest Fourier Transform in the West, version 3
URL = http://fftw.org/
DOCURL = http://www.fftw.org/fftw3_doc/Installation-on-Unix.html
PACKAGEVERSION = 3.3.10
DOWNLOADURL = http://fftw.org/fftw-${PACKAGEVERSION}.tar.gz
MODE = mpi

MODE = mpi

BUILDSYSTEM = cmake

MODULES = 

CMAKEFLAGS = \
  -D ENABLE_OMP=ON \
  -D BUILD_TESTS=OFF
CONFIGUREFLAGS = \
    --with-pic \
    --enable-shared \
    --enable-openmp \
    --enable-threads \
    --disable-dependency-tracking \
    --enable-mpi

## extra flags for intel
COMPILER==intel CONFIGUREFLAGS += --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 

PKGCONFIGLIB = pkgconfig
PREFIXPATHSET = 1
