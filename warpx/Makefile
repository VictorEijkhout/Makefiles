################################################################
####
#### Makefile for WarpX
####
################################################################

PACKAGE = warpx
PACKAGEVERSION = git
MODE = mpi
MODULES =
#blaspp lapackpp

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

GITREPO= https://github.com/ECP-WarpX/WarpX.git
include ${MAKEINCLUDES}/Make.git
