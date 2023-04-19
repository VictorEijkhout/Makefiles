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
PACKAGEVERSION = 1.14
MODE = mpi
MODULES = zlib

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
  -D HDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON \
  -D ZLIB_DIR:PATH=${TACC_ZLIB_DIR} \
  -D ZLIB_INCLUDE_DIR:PATH=${TACC_ZLIB_INC} \
  -D ZLIB_LIBRARY:PATH=${TACC_ZLIB_LIB}/libz.so \
  \
  -D HDF_ENABLE_PARALLEL:BOOL=ON \
  -D CMAKE_C_FLAGS=--enable-parallel \
  -D CMAKE_CXX_FLAGS=--enable-parallel

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
PKGCONFIGSET = lib/pkgconfig
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${PACKAGEVERSION}/hdf5-${PACKAGEVERSION}.0/src/hdf5-${PACKAGEVERSION}.0.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

