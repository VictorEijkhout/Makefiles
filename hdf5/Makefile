################################################################
####
#### Makefile for Hdf5 installation
####
#### https://hdf5.geodynamics.org/
#### https://github.com/geodynamics/hdf5
####
################################################################

PACKAGE = HDF5
PACKAGEVERSION = 1.14
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
  -D HDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON \
  -D ZLIB_INCLUDE_DIR:PATH=${TACC_ZLIB_INC} \
  -D ZLIB_LIBRARY:FILEPATH=${TACC_ZLIB_LIB} 

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
PKGCONFIGSET = lib/pkgconfig
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${PACKAGEVERSION}/hdf5-${PACKAGEVERSION}.0/src/hdf5-${PACKAGEVERSION}.0.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

