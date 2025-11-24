################################################################
####
#### Makefile for scalapack
####
################################################################

PACKAGE = SCALAPACK
ABOUT  = Scalable parallel lapack
PACKAGEVERSION = git
# 2.2.2
# release can not be compiled, waitinf for new:
# https://github.com/Reference-ScaLAPACK/scalapack.git
URL = https://github.com/Reference-ScaLAPACK/scalapack
MODE = mpi

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = 
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
NOINC = 1
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build
# requires final editing of scalapack.pc
# to have an empty "Requires:" line

GITREPO = https://github.com/Reference-ScaLAPACK/scalapack.git
include ${MAKEINCLUDES}/Make.git
DOWNLOADURL = https://github.com/Reference-ScaLAPACK/scalapack/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
