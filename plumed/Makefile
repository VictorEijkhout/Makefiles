################################################################
####
#### Makefile for Plumed installation
####
#### even the release seems to want to submodule
####
################################################################

PACKAGE = plumed
URL = https://github.com/plumed/plumed2
PACKAGEVERSION = 2.9.3

MODE = mpi
MODULES = fftw3 gsl zlib \
    mkl

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

include ${MAKEINCLUDES}/Make.configure
# most of these enable flags are yes by default;
# configure fails gracefully if the package is not found
# blas/lapack ticket:
# https://github.com/plumed/plumed2/issues/1279
CONFIGUREFLAGS = \
    --enable-external-blas=yes \
    --enable-external-lapack=yes \
    --enable-fftw=yes \
    --enable-gsl=yes \
    --enable-python=yes \
    --enable-zlib=yes
CFLAGS   = -I${TACC_FFTW3_INC} -I${TACC_GSL_INC} -I${TACC_ZLIB_INC} -I${TACC_MKL_INC}
CXXFLAGS = -I${TACC_FFTW3_INC} -I${TACC_GSL_INC} -I${TACC_ZLIB_INC} -I${TACC_MKL_INC}
CPPFLAGS = -I${TACC_FFTW3_INC} -I${TACC_GSL_INC} -I${TACC_ZLIB_INC} -I${TACC_MKL_INC}
CONFIGUREENV = PYTHON_BIN=$$( which python3 )
LIBS = \
    -L${TACC_FFTW3_LIB} -lfftw3 -L${TACC_GSL_LIB} -lgsl -L${TACC_ZLIB_LIB} -lz \
    -L${TACC_MKL_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread 

include ${MAKEINCLUDES}/Make.install
HASBIN = 1

TGZURL = https://github.com/plumed/plumed2/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

