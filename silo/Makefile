################################################################
####
#### Makefile for Silo installation
####
################################################################

BUILD = 4
# or BUILD=git

PACKAGE = SILO
URL = http://software.llnl.gov/Silo/
ifeq "${BUILD}" "git"
  PACKAGEVERSION = git
else
  PACKAGEVERSION = 4.11.1
endif
## WEIRD main branch has an error that is fixed in 4.11
## https://github.com/LLNL/Silo/issues/375
MODE = seq
MODULES = hdf5

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D SILO_ENABLE_HDF5=ON \
    -D SILO_ENABLE_FORTRAN=ON
ifeq "${TACC_FAMILY_COMPILER}" "intel"
    CMAKEFLAGS += -D CMAKE_C_FLAGS=-Wno-int-conversion
endif

HASBIN=1
CMAKEPREFIXPATHSET = 1
#DEPENDSON = hdf5
## silo executables need exact .so version
DEPENDSONCURRENT = hdf5
ifeq "${BUILD}" "git"
  include ${MAKEINCLUDES}/Make.cmake
  include ${MAKEINCLUDES}/Make.cbuild
else
  include ${MAKEINCLUDES}/Make.configure
  include ${MAKEINCLUDES}/Make.install
endif

info ::
	@echo "make default_install : configure build"
default_install : configure build

##
## download
##
TXZURL = https://github.com/LLNL/Silo/releases/download/${PACKAGEVERSION}/silo-${PACKAGEVERSION}.tar.xz
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
