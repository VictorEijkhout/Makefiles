################################################################
####
#### Makefile for Netcdf installation
#### C basis version
####
################################################################

PACKAGE = NETCDFC
URL = https://github.com/Unidata/netcdf-c
PACKAGEVERSION = 4.9.2
MODE = mpi
MODULES = phdf5

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install :
	@make configure build JCOUNT=${JCOUNT}

DOWNLOADURL = https://github.com/Unidata/netcdf-c/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
# https://downloads.unidata.ucar.edu/netcdf-c/${PACKAGEVERSION}/netcdf-c-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

