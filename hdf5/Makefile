################################################################
####
#### Makefile for Hdf5 installation
####
#### note: also comes in C++ and Fortran download
#### 
#### see release_docs/INSTALL_CMake.txt
#### 
#### nvfortran real2 support missing
#### https://hdfgroup.atlassian.net/servicedesk/customer/portal/6/HELP-2790?created=true
#### 
#### cmake DEPENDS keyword
#### https://github.com/HDFGroup/hdf5/issues/5682
#### 
################################################################

PACKAGE = HDF5
ABOUT =  Hierarchical Data Format, version 5
URL = https://github.com/HDFGroup/hdf5

# there is a micro version but we ignore that for now.
PACKAGEVERSION = 2.0.0
PACKAGESHORTVERSION = 2.0
# 1.14.6
# 1.14
## ? PACKAGEVERSIONDIR = 2_0_0

MODE = mpi
MODULES = zlib

# there is something wrong with the szip installation:
# CMake Error at /Users/eijkhout/Installation/szip/build-2.1.1-macbookair-gcc14/szip-config.cmake:61 (include):
#   include could not find requested file:

#     /Users/eijkhout/Installation/share/cmake/szip-targets.cmake


include ${MAKEINCLUDES}/Make.info

# parallel & c++ incompatible?!
# set the ALLOW_UNSUPPORTED option to fix
# in version 2 this will become HDF5_ALLOW_UNSUPPORTED
# 
HDF5PARALLEL = ON
HDFFORTRAN = ON

#### test failure with NVidia compiler:
#### https://hdfgroup.atlassian.net/servicedesk/customer/portal/6/HELP-2539
TESTING = ON

THREADSAFE=ON
CMAKEFLAGS = \
  -D HDF5_ENABLE_PLUGIN_SUPPORT=OFF \
  -D HDF5_ENABLE_THREADSAFE=${THREADSAFE} \
  -D BUILD_TESTING:BOOL=${TESTING} \
  \
  -D HDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON \
  -D ZLIB_DIR:PATH=${TACC_ZLIB_DIR} \
  -D ZLIB_INCLUDE_DIR:PATH=${TACC_ZLIB_INC} \
  -D ZLIB_LIBRARY:PATH=${TACC_ZLIB_LIB}/libz.so \
  \
  -D ALLOW_UNSUPPORTED=ON \
  -D HDF5_ALLOW_UNSUPPORTED=ON \
  -D HDF5_BUILD_FORTRAN:BOOL=${HDFFORTRAN} \
  -D HDF5_BUILD_CPP_LIB:BOOL=ON \
  -D HDF5_ENABLE_PARALLEL:BOOL=${HDF5PARALLEL}
ifeq "${TACC_FAMILY_COMPILER}" "gcc"
  CMAKEFLAGS += -DCMAKE_EXE_LINKER_FLAGS='-lpthread -lgfortran'
else
  CMAKEFLAGS += -DCMAKE_EXE_LINKER_FLAGS='-lpthread'
endif

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
MODULENOTES = Thread-safe mode enabled
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
	@if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then \
	    export I_MPI_FCFLAGS="-I${TACC_IMPI_INC}/mpi/gfortran/11.1.0 -I${TACC_IMPI_INC}/gfortran/11.1.0" \
	 ; fi \
	 && make --no-print-directory configure build JCOUNT=${JCOUNT} \
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

# https://support.hdfgroup.org/releases/hdf5/v2_0/v2_0_0/downloads/index.html

# version 2 semantic versioning at last. and they still mess it up
## TGZURL = https://support.hdfgroup.org/releases/hdf5/v2_0/v2_0_0/downloads/hdf5-${PACKAGEVERSION}.tar.gz

# my own attempt at predictable naming
TGZURL = https://github.com/VictorEijkhout/hdf5-with-predictable-naming/archive/refs/tags/${PACKAGEVERSION}.tar.gz

include ${MAKEINCLUDES}/Make.download

##
## git clone
##
GITREPO = https://github.com/HDFGroup/hdf5.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean
