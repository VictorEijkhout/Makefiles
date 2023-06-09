################################################################
####
#### Makefile for Hdf5 installation
####
#### note: also comes in C++ and Fortran download
#### https://portal.hdfgroup.org/display/HDF5/Parallel+HDF5
####
################################################################

PACKAGE = HDF5
# there is a micro version but we ignore that for now.
PACKAGESHORTVERSION = 1.14
PACKAGEVERSION = 1.14.0
# 1.14.1 is out but uses weird-ass naming of the tgz
MODE = mpi
MODULES = zlib

include ${MAKEINCLUDES}/Make.info

HDF5PARALLEL = ON
VARIANT = par
CMAKEFLAGS = \
  -D HDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON \
  -D ZLIB_DIR:PATH=${TACC_ZLIB_DIR} \
  -D ZLIB_INCLUDE_DIR:PATH=${TACC_ZLIB_INC} \
  -D ZLIB_LIBRARY:PATH=${TACC_ZLIB_LIB}/libz.so \
  \
  -D HDF5_ENABLE_PARALLEL:BOOL=${HDF5PARALLEL}

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
PKGCONFIG = lib/pkgconfig
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make par : parallel installation"
	@echo "make seq : sequential installation"
.PHONY: seq par
seq :
	make configure build VARIANT=seq HDF5PARALLEL=OFF MODE=seq
par :
	make configure build VARIANT=par HDF5PARALLEL=ON

TGZURL = https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${PACKAGESHORTVERSION}/hdf5-${PACKAGEVERSION}/src/hdf5-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

