################################################################
####
#### Makefile for Octopus installation
####
################################################################

PACKAGE = Octopus
URL = https://gitlab.com/octopus-code/octopus \
    https://www.octopus-code.org/documentation/15/manual/installation/cmake/
PACKAGEVERSION = 16.0
GITTAG = 16.0

MODE = mpi
MODULES = gsl mkl 
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
include ${MAKEINCLUDES}/Make.cmake
CMAKEFLAGS = \
    -D CMAKE_DISABLE_FIND_PACKAGE_Libxc=On \
    -D CMAKE_DISABLE_FIND_PACKAGE_METIS=On \
    -D OCTOPUS_OpenMP=ON \
    -D OCTOPUS_MPI=${OCTOPUS_MPI} \
    -D OCTOPUS_MKL=ON  \
    -D OCTOPUS_ScaLAPACK=ON
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
