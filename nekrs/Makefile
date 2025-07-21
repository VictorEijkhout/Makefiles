################################################################
####
#### Makefile for NekRS installation
####
################################################################

PACKAGE = nekRS
URL = https://github.com/Nek5000/nekRS
PACKAGEVERSION = 23.0

MODE = mpi
MODULES_BASIC = 
MODULES = ${MODULES_BASIC} 

include ${MAKEINCLUDES}/Make.info

##
## Cmake 
##
USEFORTRAN = ON
CMAKEFLAGS_BASIC = 
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
	  CMAKEFLAGS="${CMAKEFLAGS_BASIC} -D NEKRS_USE_CUDA:BOOL=ON"

##
## Download
##
TGZURL = https://github.com/Nek5000/nekRS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## Clean
##
include ${MAKEINCLUDES}/Make.clean

