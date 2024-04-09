################################################################
####
#### Makefile for Hdf5 installation
####
#### note: also comes in C++ and Fortran download
#### 
#### see release_docs/INSTALL_CMake.txt
#### 
################################################################

PACKAGE = HDF5
URL = https://portal.hdfgroup.org/downloads/index.html
# there is a micro version but we ignore that for now.
PACKAGESHORTVERSION = 1.14
PACKAGEVERSION = 1.14.3
MODE = mpi
MODULES = zlib

include ${MAKEINCLUDES}/Make.info

HDF5PARALLEL = ON
ifeq "${TACC_FAMILY_COMPILER}" "gcc"
  HDFFORTRAN = OFF
else
  HDFFORTRAN = ON
endif

ifeq "${TACC_SYSTEM}" "macbookair"
  LIBEXT = dylib
else
  LIBEXT = SO
endif
CMAKEFLAGS = \
  -D BUILD_SHARED_LIBS=TRUE \
  -D HDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON \
  -D ZLIB_DIR:PATH=${TACC_ZLIB_DIR} \
  -D ZLIB_INCLUDE_DIR:PATH=${TACC_ZLIB_INC} \
  -D ZLIB_LIBRARY:PATH=${TACC_ZLIB_LIB}/libz.${LIBEXT} \
  \
  -D HDF5_BUILD_FORTRAN:BOOL=${HDFFORTRAN} \
  -D HDF5_ENABLE_PARALLEL:BOOL=${HDF5PARALLEL}

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make par : parallel installation"
	@echo "make seq : sequential installation"
.PHONY: seq par
seq :
	@make configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    HDF5PARALLEL=OFF MODE=seq MODULENAME=hdf5
par :
	@make configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    HDF5PARALLEL=ON  MODE=mpi MODULENAME=phdf5 MODULENAMEALT=hdf5
info ::
	@echo "make default_install : sequential and parallel versions"
default_install :: 
	@make --no-print-directory seq par JCOUNT=${JCOUNT} PACKAGEVERSION=${PACKAGEVERSION} \
	    seq par

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
