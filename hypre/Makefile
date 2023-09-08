################################################################
####
#### Makefile for Hypre installation
####
################################################################

PACKAGE = hypre
URL = https://github.com/hypre-space/hypre
PACKAGEVERSION = 2.29.0

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

.PHONY: small big
small :
	@make --no-print-directory configure build \
	    MODULENAME=hypre32 MODULENAMEALT=hypre
big :
	@make --no-print-directory configure build \
	    CONFIGUREFLAGS=--enable-bigint \
	    MODULENAME=hypre64 MODULENAMEALT=hypre

CONFIGURESUBDIR = src
MAKESUBDIR = src
PKGCONFIGLIB = pkgconfig

##
## configure version
##
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
##
## cmake version
##
# include ${MAKEINCLUDES}/Make.cmake
# include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/hypre-space/hypre/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

