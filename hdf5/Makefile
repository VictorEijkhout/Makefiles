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
# https://github.com/HDFGroup/hdf5
# https://support.hdfgroup.org/ftp/HDF5/releases/

# there is a micro version but we ignore that for now.
PACKAGESHORTVERSION = 1.14
PACKAGEVERSION = 1.14.6

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

CMAKEFLAGS = \
  -D HDF5_ENABLE_PLUGIN_SUPPORT=OFF \
  -D HDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON \
  -D BUILD_TESTING:BOOL=${TESTING} \
  -D ZLIB_DIR:PATH=${TACC_ZLIB_DIR} \
  -D ZLIB_INCLUDE_DIR:PATH=${TACC_ZLIB_INC} \
  -D ZLIB_LIBRARY:PATH=${TACC_ZLIB_LIB}/libz.so \
  \
  -D ALLOW_UNSUPPORTED=ON \
  -D HDF5_ALLOW_UNSUPPORTED=ON \
  -D HDF5_BUILD_FORTRAN:BOOL=${HDFFORTRAN} \
  -D HDF5_BUILD_CPP_LIB:BOOL=ON \
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
	@if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then \
	    export I_MPI_FCFLAGS="-I${TACC_IMPI_INC}/mpi/gfortran/11.1.0" \
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

# my own attempt at predictable naming
TGZURL = https://github.com/VictorEijkhout/hdf5-with-predictable-naming/archive/refs/tags/${PACKAGEVERSION}.tar.gz

# if hdf group ever gets semantic versioning straightened out
## https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${PACKAGESHORTVERSION}/hdf5-${PACKAGEVERSION}/src/hdf5-${PACKAGEVERSION}.tar.gz

include ${MAKEINCLUDES}/Make.download

##
## git clone
##
GITREPO = https://github.com/HDFGroup/hdf5.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean
