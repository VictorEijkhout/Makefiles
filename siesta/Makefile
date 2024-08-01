################################################################
####
#### Makefile for Siesta installation
####
#### https://gitlab.com/siesta-project/siesta/-/issues/397
####
################################################################

PACKAGE = siesta
URL = https://siesta-project.org/siesta/index.html
PACKAGEVERSION = 5.0.1

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

##
## cmake version
##
# https://docs.siesta-project.org/projects/siesta/en/stable/installation/build-issues.html#siesta-does-not-recognize-intel-mkl-libraries
CMAKEFLAGS = \
    -D SIESTA_WITH_FLOOK=FALSE \
    -D Fortran_FLAGS_DEBUG=-O0 \
    -D SCALAPACK_INCLUDE=${TACC_INTEL_DIR}/include/mkl_scalapack.h \
    -D SCALAPACK_LIBRARY_DIR=${TACC_INTEL_LIB} \
    -D SCALAPACK_LIBRARY='-lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_lp64'
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build


TGZURL = https://gitlab.com/siesta-project/siesta/-/archive/${PACKAGEVERSION}/siesta-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://gitlab.com/siesta-project/siesta.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git
