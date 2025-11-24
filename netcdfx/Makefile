################################################################
####
#### Makefile for Netcdf installation
#### C++ interface
####
################################################################

PACKAGE = NETCDFX
URL = https://github.com/Unidata/netcdf-cxx4
PACKAGEVERSION = 4.3.1
MODE = seq
MODULES =
# hdf5 netcdf
NONMODULES =
# phdf5 parallelnetcdf

include ${MAKEINCLUDES}/Make.info

#
# cmake discovery seems broken, hence explicit lib/inc
# the hdf5 lib variable appears in the cmakelists,
# the hdf5 inc is pure guesswork
# https://github.com/Unidata/netcdf-cxx4/issues/151
#
CMAKEFLAGS = \
    -D HDF5_C_LIBRARY_hdf5=${TACC_HDF5_LIB}/libhdf5.so \
    -D NOT_HDF5_INCLUDE_DIRECTORIES=${TACC_HDF5_INC} \
    -D netCDF_LIBRARIES=${TACC_NETCDF_LIB}/libnetcdf.so \
    -D netCDF_INCLUDE_DIR='${TACC_NETCDF_INC};${TACC_HDF5_INC}'
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
# this fixes the compile problem
# but now the h5 include files wind up in the netcdf inc
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

DOWNLOADURL = https://github.com/Unidata/netcdf-cxx4/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean
