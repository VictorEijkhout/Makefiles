################################################################
####
#### Makefile for Fftw3 installation
####
################################################################

PACKAGE = FFTW3
ABOUT = Fastest Fourier Transform in the West, version 3
URL = http://fftw.org/
DOCURL = http://www.fftw.org/fftw3_doc/Installation-on-Unix.html
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

## insert: AC_DEFUN([_AC_PROG_FC_V], [])
DEFUNPROGFC = 1
FORCERECONF = 1

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
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
	    \
	    MODE=mpi
## INSTALLVARIANT=single INSTALLEXT=single 
## MODULEVERSIONEXTRA=single
double :
	@make configure build JCOUNT=${JCOUNT} \
	    \
	    MODE=mpi
## INSTALLVARIANT=double 

DOWNLOADURL = http://fftw.org/fftw-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

