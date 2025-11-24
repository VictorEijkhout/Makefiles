################################################################
####
#### Makefile for Octopus installation
#### autotools install
####
################################################################

PACKAGE = Octopus
URL = https://gitlab.com/octopus-code/octopus \
    https://www.octopus-code.org/documentation/15/manual/installation/cmake/
PACKAGEVERSION = 15.1

MODE = mpi
MODULES = gsl fftw3 libxc mkl

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

include ${MAKEINCLUDES}/Make.configure
BEFORECONFIGURECMDS = autoupdate
CONFIGUREFLAGS = \
    -enable-mpi --enable-openmp \
    --with-blas="-L${TACC_MKL_LIB} -lmkl_rt -lpthread -lm -ldl" \
    --with-libxc-prefix=${TACC_LIBXC_DIR} \
    --with-fftw-prefix=${TACC_FFTW3_DIR} \
    --with-blacs="-L${TACC_MKL_LIB} -lmkl_rt -lpthread -lm -ldl" \
    --with-scalapack="-L${TACC_MKL_LIB} -lmkl_rt -lpthread -lm -ldl" \
    --with-gsl-prefix=${TACC_GSL_DIR}
include ${MAKEINCLUDES}/Make.install
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

GITREPO = https://gitlab.com/octopus-code/octopus.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git

DOWNLOADURL = https://gitlab.com/octopus-code/octopus/-/archive/${PACKAGEVERSION}/octopus-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

