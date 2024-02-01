################################################################
####
#### Makefile for Netcdf installation
#### Fortran interface
####
################################################################

PACKAGE = NETCDFF
URL = https://github.com/Unidata/netcdf-fortran \
    https://docs.unidata.ucar.edu/netcdf-fortran/current \
    https://docs.unidata.ucar.edu/netcdf-c/current/building_netcdf_fortran.html
PACKAGEVERSION = 4.6.1
MODE = mpi
MODULES = phdf5 netcdf

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = -DNCDIR=${TACC_NETCDF_DIR}
include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install :
	@make configure build JCOUNT=${JCOUNT}

TGZURL = https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

