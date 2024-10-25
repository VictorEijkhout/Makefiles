################################################################
####
#### Makefile for Charm++ installation
####
################################################################

PACKAGE = charmpp
URL = https://charmplusplus.org/
PACKAGEVERSION = 8.0.1

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

CMAKEFLAGS = \
    -D NETWORK=mpi \
    -D CMAKE_POSITION_INDEPENDENT_CODE=ON
BUILDSTATICLIBS = 1

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build


TGZURL = https://github.com/charmplusplus/charm/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

