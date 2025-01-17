################################################################
####
#### Makefile for scalapack
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
# requires final editing of scalapack.pc
# to have an empty "Requires:" line

GITREPO = https://github.com/Reference-ScaLAPACK/scalapack.git
TGZURL = https://github.com/Reference-ScaLAPACK/scalapack/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
