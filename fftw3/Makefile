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
CONFIGUREFLAGS = \
    --with-pic \
    --enable-shared \
    --enable-openmp \
    --enable-threads \
    --disable-dependency-tracking \
    --enable-mpi
ifeq "${TACC_FAMILY_COMPILER}" "intel"
CONFIGUREFLAGS += \
    --enable-sse2 \
    --enable-avx \
    --enable-avx2 \
    --enable-avx512 
endif

PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
default_install : configure build

TGZURL = http://fftw.org/fftw-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

