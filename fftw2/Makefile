################################################################
####
#### Makefile for Fftw3 installation
####
################################################################

PACKAGE = FFTW2
ABOUT = Fastest Fourier Transform in the West, version 2
URL = http://fftw.org/
PACKAGEVERSION = 2.1.5
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
CONFIGUREFLAGS_STD = \
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
CONFIGUREFLAGS = ${CONFIGUREFLAGS_STD}

## PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
info ::
	@echo "make single, double"
	@echo "make default_install : double precision version"
default_install : double
.PHONY: single double
single :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    CONFIGUREFLAGS="${CONFIGUREFLAGS_STD} --enable-float" \
	    INSTALLVARIANT=single INSTALLEXT=single
double :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    INSTALLVARIANT=double 

TGZURL = http://fftw.org/fftw-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

## config.guess can not guess NVIDIA
## http://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.guess;hb=HEAD
## http://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.sub;hb=HEAD
