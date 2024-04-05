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
MODE = seq
MODULES = hdf5 netcdf
NONMODULES = phdf5 parallelnetcdf

include ${MAKEINCLUDES}/Make.info

#
# 1. NCDIR not actually used though it appears in the docs
# 2. cmake discovery seems broken, hence explicit lib/inc
#
CMAKEFLAGS = -DNCDIR=${TACC_NETCDF_DIR} \
    -DnetCDF_LIBRARIES=${TACC_NETCDF_LIB}/libnetcdf.so \
    -DnetCDF_INCLUDE_DIR=${TACC_NETCDF_INC}
include ${MAKEINCLUDES}/Make.cmake
##
## Install into netcdf installation
## 
INSTALLPATH = ${TACC_NETCDF_DIR}
NOMODULE = 1
#HASBIN=1
#PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: seq par default_install
default_install : seq
info ::
	@echo "make seq par, default_install=seq"
seq : configure build
par :
	@make configure build JCOUNT=${JCOUNT} \
	  MODULES="phdf5 parallelnetcdf" \
	  NONMODULES="hdf5 netcdf" \
	  CMAKEFLAGS="\
	    -DnetCDF_LIBRARIES=${TACC_PARALLELNETCDF_LIB}/libpnetcdf.so \
	    -DnetCDF_INCLUDE_DIR=${TACC_PARALLELNETCDF_INC} \
	    " \
	  INSTALLPATH=${TACC_PARALLELNETCDF_DIR}

TGZURL = https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

