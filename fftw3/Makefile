################################################################
####
#### Makefile for Fftw3 installation
####
################################################################

PACKAGE = FFTW3
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 3.3.10
MAKEINCLUDES = ${STOCKYARD}/makefiles
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
INSTALL_USE_DESTDIR  = true
## make DESTDIR=$RPM_BUILD_ROOT install

TGZURL = http://fftw.org/fftw-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
