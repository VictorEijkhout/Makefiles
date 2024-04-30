################################################################
####
#### Makefile for Fftw3 installation
####
################################################################

PACKAGE = FFTW3
URL = http://fftw.org/ \
    http://www.fftw.org/fftw3_doc/Installation-on-Unix.html
PACKAGEVERSION = 3.3.10
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

PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.install
.PHONY: default_install
info ::
	@echo "================ fftw rules"
	@echo "make single double"
	@echo "make default_install : double precision version"
default_install : double
.PHONY: single double
single :
	@make configure build JCOUNT=${JCOUNT} \
	    CONFIGUREFLAGS="${CONFIGUREFLAGS_STD} --enable-float" \
	    INSTALLVARIANT=single INSTALLEXT=single \
	    MODE=mpi
double :
	@make configure build JCOUNT=${JCOUNT} \
	    INSTALLVARIANT=double \
	    MODE=mpi

TGZURL = http://fftw.org/fftw-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

