################################################################
####
#### Makefile for Octopus installation
####
################################################################

PACKAGE = Octopus
URL = https://gitlab.com/octopus-code/octopus
DOCURL = https://www.octopus-code.org/documentation/15/manual/installation/cmake/
PACKAGEVERSION = 16.0
GITTAG = 16.0

MODE = mpi
MODULES = gsl scalapack \
    mkl 
ifeq "${TACC_FAMILY_COMPILER}" "gcc"
  MODULES += ninja
endif
#  fftw3 not needed because MKL

# Note: metis is included https://www.octopus-code.org/documentation/15/faq/
# Issue: https://gitlab.com/octopus-code/octopus/-/issues/1157

# Git clone fails on libxc:
# https://gitlab.com/octopus-code/octopus/-/issues/1158
# seems resolved by disable finding Libxc

include ${MAKEINCLUDES}/Make.info

##
## cmake install (there is an autoconf one)
## https://www.octopus-code.org/documentation/15/manual/installation/cmake/
##
OCTOPUS_MPI = ON
ifeq "${TACC_SYSTEM}" "vista"
  # cmake 4.1 should be able to find nvpl
  USEMKL=OFF
  DEFINE_LAPACK="-DBLA_VENDOR=NVPL"
else
  USEMKL=ON
  DEFINE_LAPACK='-D OCTOPUS_MKL=ON'
endif
include ${MAKEINCLUDES}/Make.cmake
CMAKEUSENINJA = 1
CMAKEFLAGS = \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -D CMAKE_DISABLE_FIND_PACKAGE_Libxc=On \
    -D CMAKE_DISABLE_FIND_PACKAGE_METIS=On \
    -D OCTOPUS_OpenMP=ON \
    -D OCTOPUS_MPI=${OCTOPUS_MPI} \
    ${DEFINE_LAPACK} \
    -D OCTOPUS_ScaLAPACK=ON \
    -D BUFFER=ON
## Note: there are also "--preset stufstufstuf" flags for cmake. See ticket 1158.

include ${MAKEINCLUDES}/Make.cbuild
HASBIN = 1
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build
DEPENDSON = gsl

GITREPO = https://gitlab.com/octopus-code/octopus.git
SUBMODULE = 1
# GITTAG with version number listed up top
include ${MAKEINCLUDES}/Make.git

#TGZURL = https://octopus-code.org/download/${PACKAGEVERSION}/octopus-${PACKAGEVERSION}.tgz
TGZURL = https://gitlab.com/octopus-code/octopus/-/archive/${PACKAGEVERSION}/octopus-${PACKAGEVERSION}.tar.gz
#        https://gitlab.com/octopus-code/octopus/-/archive/16.0/octopus-16.0.tar.gz
# this is the downloaded tar ball, which can not be used 
# TXZURL = https://octopus-code.org/download/${PACKAGEVERSION}/octopus-${PACKAGEVERSION}.tar.xz
include ${MAKEINCLUDES}/Make.download
