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
URL = https://portal.hdfgroup.org/downloads/index.html \
    https://support.hdfgroup.org/ftp/HDF5/releases/
# there is a micro version but we ignore that for now.
PACKAGESHORTVERSION = 1.14
PACKAGEVERSION = 1.14.4
MODE = mpi
MODULES = zlib

# there is something wrong with the szip installation:
# CMake Error at /Users/eijkhout/Installation/szip/build-2.1.1-macbookair-gcc14/szip-config.cmake:61 (include):
#   include could not find requested file:

#     /Users/eijkhout/Installation/share/cmake/szip-targets.cmake


include ${MAKEINCLUDES}/Make.info

HDF5PARALLEL = ON
HDFFORTRAN = ON

ifeq "${TACC_SYSTEM}" "macbookair"
  LIBEXT = dylib
else
  LIBEXT = so
endif

#### test failure with NVidia compiler:
#### https://hdfgroup.atlassian.net/servicedesk/customer/portal/6/HELP-2539
TESTING = ON

CMAKEFLAGS = \
  -D HDF5_ENABLE_PLUGIN_SUPPORT=OFF \
  -D HDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON \
  -D BUILD_TESTING:BOOL=${TESTING} \
  -D ZLIB_DIR:PATH=${TACC_ZLIB_DIR} \
  -D ZLIB_INCLUDE_DIR:PATH=${TACC_ZLIB_INC} \
  -D ZLIB_LIBRARY:PATH=${TACC_ZLIB_LIB}/libz.${LIBEXT} \
  \
  -D HDF5_BUILD_FORTRAN:BOOL=${HDFFORTRAN} \
  -D HDF5_ENABLE_PARALLEL:BOOL=${HDF5PARALLEL}

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
FAMILY = hf5
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "================ hdf5 rules"
	@echo "make par : parallel installation"
	@echo "make seq : sequential installation"
.PHONY: seq par
seq :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    HDF5PARALLEL=OFF MODE=seq HDFFORTRAN=${HDFFORTRAN} \
	    MODULENAME=hdf5
par :
	@make --no-print-directory configure build JCOUNT=${JCOUNT} \
	    PACKAGEVERSION=${PACKAGEVERSION} \
	    HDF5PARALLEL=ON  MODE=mpi HDFFORTRAN=${HDFFORTRAN} \
	    MODULENAME=phdf5 MODULENAMEALT=hdf5
info ::
	@echo "make default_install : sequential and parallel versions"
default_install :: 
	@make --no-print-directory seq par JCOUNT=${JCOUNT} PACKAGEVERSION=${PACKAGEVERSION} \
	    seq par

##
## download
##
TGZURL = https://github.com/VictorEijkhout/hdf5-with-predictable-naming/archive/refs/tags/${PACKAGEVERSION}.tar.gz
## https://github.com/VictorEijkhout/hdf5-with-predictable-naming/releases/download/${PACKAGEVERSION}/hdf5-${PACKAGEVERSION}.tar.gz
## https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${PACKAGESHORTVERSION}/hdf5-${PACKAGEVERSION}/src/hdf5-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## git clone
##
GITREPO = https://github.com/HDFGroup/hdf5.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean
