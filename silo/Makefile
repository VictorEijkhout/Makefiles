################################################################
####
#### Makefile for Silo installation
####
#### note: also comes in C++ and Fortran download
#### 
#### see release_docs/INSTALL_CMake.txt
#### 
################################################################

PACKAGE = SILO
URL = http://software.llnl.gov/Silo/
PACKAGEVERSION = git
#4.11.1
MODE = seq
MODULES = hdf5

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D SILO_ENABLE_HDF5=ON \
    -D SILO_ENABLE_FORTRAN=ON
HASBIN=1
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
# include ${MAKEINCLUDES}/Make.configure
# include ${MAKEINCLUDES}/Make.install
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
include ${MAKEINCLUDES}/Make.clean
