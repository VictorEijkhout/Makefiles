################################################################
####
#### Makefile for Adios2 installation
####
################################################################

PACKAGE = adios2
URL = https://github.com/ornladios/ADIOS2
PACKAGEVERSION = 2.10.2

MODE = mpi
MODULES_BASIC = phdf5 
MODULES = ${MODULES_BASIC} 
## can use hdf5, kokkos, zlib

include ${MAKEINCLUDES}/Make.info

## PKGCONFIGLIB = pkgconfig

##
## cmake 
##
CMAKEFLAGS_BASIC = -D CMAKE_BUILD_TYPE=Release \
    -D ADIOS2_USE_Fortran=OFF \
    -D ADIOS2_USE_AWSSDK=OFF
CMAKEFLAGS = ${CMAKEFLAGS_BASIC}
HASBIN = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : cpu; also: gpu"
.PHONY: default_install cpu gpu
default_install : cpu
cpu : configure build
gpu :
	@make --no-print-directory configure build \
	  MODULES="${MODULES_BASIC} cuda" INSTALLVARIANT=cuda MODULEVERSIONEXTRA=cuda \
	  CMAKEFLAGS="${CMAKEFLAGS_BASIC} -D ADIOS2_USE_CUDA:BOOL=ON"

TGZURL = https://github.com/ornladios/ADIOS2/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = git@github.com:ornladios/ADIOS2.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

