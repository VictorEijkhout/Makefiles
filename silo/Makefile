################################################################
####
#### Makefile for Silo installation
####
################################################################

PACKAGE = SILO
ABOUT = File-based, scientific data exchange and software interoperability
URL = http://software.llnl.gov/Silo/
PACKAGEVERSION = 4.12.0
## PACKAGEVERSION_DOWNLOAD = 4.12.0-pre1

## WEIRD main branch has an error that is fixed in 4.11
## https://github.com/LLNL/Silo/issues/375
## Intel 24 problem resolved?
## https://github.com/LLNL/Silo/issues/400
MODE = seq
MODULES = hdf5

include ${MAKEINCLUDES}/Make.info

HASBIN=1
CMAKEPREFIXPATHSET = 1
DEPENDSONCURRENT = hdf5

## starting with 4.12 only cmake build
CMAKEFLAGS = \
      -D SILO_ENABLE_HDF5=ON \
      -D SILO_ENABLE_FORTRAN=ON

ifeq "${TACC_FAMILY_COMPILER}" "intel"
      CMAKEFLAGS += -D CMAKE_C_FLAGS=-Wno-int-conversion
endif

##
## This is needed on Frontera, not elsewhere
##
ifeq "${TACC_FAMILY_COMPILER}" "gcc"
  CMAKEFLAGS += -DCMAKE_EXE_LINKER_FLAGS='-lpthread -lgfortran'
else
  CMAKEFLAGS += -DCMAKE_EXE_LINKER_FLAGS='-lpthread'
endif

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
default_install : configure build

##
## download
##
DOWNLOADURL = https://github.com/LLNL/Silo/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

##
## git clone
##
GITREPO = https://github.com/LLNL/Silo.git
include ${MAKEINCLUDES}/Make.git

##
## Clean
##
include ${MAKEINCLUDES}/Make.clean
