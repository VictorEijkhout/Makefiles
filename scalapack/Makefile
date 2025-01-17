################################################################
####
#### Makefile for scalapack
####
#### Note: Conquest only supports version 4 and 5, but version 6 is out
####
################################################################

PACKAGE = SCALAPACK
PACKAGEVERSION = 2.2.2
URL = https://github.com/Reference-ScaLAPACK/scalapack
MODE = mpi

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = 
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

GITREPO = https://github.com/Reference-ScaLAPACK/scalapack.git
TGZURL = https://github.com/Reference-ScaLAPACK/scalapack/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
