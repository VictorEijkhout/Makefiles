################################################################
####
#### Makefile for Adios2 installation
####
################################################################

PACKAGE = adios2
URL = https://github.com/ornladios/ADIOS2
PACKAGEVERSION = 2.9.1

MODE = mpi
MODULES = phdf5 
## can use hdf5, kokkos, zlib

include ${MAKEINCLUDES}/Make.info

## PKGCONFIGLIB = pkgconfig

##
## cmake 
##
CMAKEFLAGS = -D CMAKE_BUILD_TYPE=Release \
    -D ADIOS2_USE_Fortran=OFF \
    -D ADIOS2_USE_AWSSDK=OFF
HASBIN = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/ornladios/ADIOS2/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = git@github.com:ornladios/ADIOS2.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

