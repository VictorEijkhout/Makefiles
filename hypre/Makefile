################################################################
####
#### Makefile for Hypre installation
####
################################################################

PACKAGE = hypre
URL = https://github.com/hypre-space/hypre \
    https://hypre.readthedocs.io/en/latest/ch-misc.html#cmake-instructions
PACKAGEVERSION = 2.31.0

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make small big ( 32/64 bit int versions)"
.PHONY: small big
small :
	@make --no-print-directory configure build public JCOUNT=${JCOUNT} \
	    CMAKEFLAGS="-D HYPRE_ENABLE_SHARED=ON" \
	    INSTALLVARIANT=i32 
big :
	@make --no-print-directory configure build public JCOUNT=${JCOUNT} \
	    CMAKEFLAGS="-D HYPRE_ENABLE_SHARED=ON -D HYPRE_ENABLE_BIGINT=ON" \
	    CONFIGUREFLAGS=--enable-bigint \
	    INSTALLVARIANT=i64 MODULEVERSIONEXTRA=i64
default_install : small big

##
## configure version
##
CONFIGURESUBDIR = src
MAKESUBDIR = src
#include ${MAKEINCLUDES}/Make.configure
#include ${MAKEINCLUDES}/Make.install

##
## cmake version
##
CMAKESUBDIR = src
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
CPTOINSTALLDIR = src/examples

TGZURL = https://github.com/hypre-space/hypre/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

