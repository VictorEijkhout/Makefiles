################################################################
####
#### Makefile for Octopus installation
####
################################################################

PACKAGE = Octopus
URL = https://gitlab.com/octopus-code/octopus \
    https://www.octopus-code.org/documentation/15/manual/installation/cmake/
PACKAGEVERSION = 15.1

MODE = mpi
MODULES = gsl fftw3 libxc mkl
# Note: metis is included https://www.octopus-code.org/documentation/15/faq/
# Issue: https://gitlab.com/octopus-code/octopus/-/issues/1157

# Git clone fails on libxc:
# https://gitlab.com/octopus-code/octopus/-/issues/1158

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

##
## cmake install (there is an autoconf one)
## https://www.octopus-code.org/documentation/15/manual/installation/cmake/
##
include ${MAKEINCLUDES}/Make.cmake
CMAKEFLAGS = \
    -D OCTOPUS_MPI=ON 
#    -D OCTOPUS_METIS=ON --debug-find-pkg=metis
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

GITREPO = https://gitlab.com/octopus-code/octopus.git
include ${MAKEINCLUDES}/Make.git
TXZURL = https://octopus-code.org/download/15.1/octopus-${PACKAGEVERSION}.tar.xz
include ${MAKEINCLUDES}/Make.download

