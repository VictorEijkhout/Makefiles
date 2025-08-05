################################################################
####
#### Makefile for Adios2 installation
####
#### uint_ types:
#### https://github.com/ornladios/ADIOS2/issues/4534
####
################################################################

PACKAGE = adios2
ABOUT = The Adaptable Input/Output System version 2
URL = https://github.com/ornladios/ADIOS2
PACKAGEVERSION = 2.10.2

MODE = mpi
MODULES_BASIC = phdf5
MODULES = ${MODULES_BASIC} 
## can use hdf5, kokkos, zlib

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
# https://adios2.readthedocs.io/en/latest/setting_up/setting_up.html#cmake-options
USEFORTRAN = ON
# --debug-find to see cmake in action
CMAKEFLAGS_BASIC =  \
    -D ADIOS2_USE_Fortran=${USEFORTRAN} \
    -D ADIOS2_USE_AWSSDK=OFF \
    -D CMAKE_EXE_LINKER_FLAGS=-lfabric
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
BOOTSTRAP = sed -i -e /iomanip/a\"\#include <cstdint>\" thirdparty/yaml-cpp/yaml-cpp/src/emitterutils.cpp
include ${MAKEINCLUDES}/Make.download
GITREPO = git@github.com:ornladios/ADIOS2.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

