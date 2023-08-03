################################################################
####
#### Makefile for Hdf5 installation
####
#### note: also comes in C++ and Fortran download
#### 
################################################################

PACKAGE = HDF5
URL = https://portal.hdfgroup.org/display/HDF5/Parallel+HDF5
# there is a micro version but we ignore that for now.
PACKAGESHORTVERSION = 1.14
PACKAGEVERSION = 1.14.0
# 1.14.1 is out but uses weird-ass naming of the tgz
MODE = mpi
MODULES = zlib

include ${MAKEINCLUDES}/Make.info

HDF5PARALLEL = ON
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
	make configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    HDF5PARALLEL=OFF MODE=seq MODULENAME=hdf5
par :
	make configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    HDF5PARALLEL=ON  MODE=mpi MODULENAME=phdf5 MODULENAMEALT=hdf5

##
## download
##
TGZURL = https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${PACKAGESHORTVERSION}/hdf5-${PACKAGEVERSION}/src/hdf5-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## git clone
##
GITREPO = https://github.com/HDFGroup/hdf5.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean

