################################################################
####
#### Makefile for Fftw3 installation
####
################################################################

PACKAGE = FFTW3
PACKAGEVERSION = 3.3.10
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
CONFIGUREFLAGS =  --with-pic \
            --enable-shared \
            --enable-openmp \
            --enable-threads \
            --disable-dependency-tracking \
            --enable-mpi \
            --enable-sse2 \
            --enable-avx \
            --enable-avx2 \
            --enable-avx512 

include ${MAKEINCLUDES}/Make.install

TGZURL = http://fftw.org/fftw-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

