################################################################
####
#### Makefile for Octopus installation
####
################################################################

PACKAGE = Octopus
URL = https://gitlab.com/octopus-code/octopus \
    https://www.octopus-code.org/documentation/15/manual/installation/cmake/
PACKAGEVERSION = 15.1
GITTAG = 15.1

MODE = mpi
MODULES = gsl fftw3 libxc mkl
# Note: metis is included https://www.octopus-code.org/documentation/15/faq/
# Issue: https://gitlab.com/octopus-code/octopus/-/issues/1157

# Git clone fails on libxc:
# https://gitlab.com/octopus-code/octopus/-/issues/1158
# seems resolved by setting Libxc_ROOT

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

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
    -D OCTOPUS_MPI=${OCTOPUS_MPI}
## Note: there are also "--preset stufstufstuf" flags for cmake. See ticket 1158.

include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

GITREPO = https://gitlab.com/octopus-code/octopus.git
SUBMODULE = 1
# GITTAG with version number listed up top
include ${MAKEINCLUDES}/Make.git

TGZURL = https://octopus-code.org/download/${PACKAGEVERSION}/octopus-${PACKAGEVERSION}.tgz
# this is the downloaded tar ball, which can not be used 
# TXZURL = https://octopus-code.org/download/${PACKAGEVERSION}/octopus-${PACKAGEVERSION}.tar.xz
include ${MAKEINCLUDES}/Make.download
