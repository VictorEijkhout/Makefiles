################################################################
####
#### Makefile for Neuralfortran installation
####
#### just testing out for user:
#### https://consult.tacc.utexas.edu/Ticket/Display.html?id=99966
####
################################################################

PACKAGE = neuralfortran
URL = https://github.com/modern-fortran/neural-fortran
PACKAGEVERSION = 0.17.0

MODE = mpi
MODULES = hdf5

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

##
## configure version
##
CMAKEFLAGS = -DSERIAL=1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install seq par
seq : configure build
default_install : seq

TGZURL = https://github.com/modern-fortran/neural-fortran/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

