################################################################
####
#### Makefile for Octopus installation
####
################################################################

PACKAGE = Octopus
URL = www.octopus-code.org
PACKAGEVERSION = 15.1

MODE = mpi
MODULES = mkl

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build


TXZURL = https://octopus-code.org/download/15.1/octopus-${PACKAGEVERSION}.tar.xz
include ${MAKEINCLUDES}/Make.download

