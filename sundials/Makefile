################################################################
####
#### Makefile for Sundials installation
####
#### https://github.com/geodynamics/sundials
####
################################################################

PACKAGE = SUNDIALS
URL = https://computing.llnl.gov/projects/sundials
PACKAGEVERSION = 7.3.0
MODE = mpi

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS =
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
info ::
	@echo "make default_install : plain configure build"
default_install : configure build

HASBIN = 1

TGZURL = https://github.com/LLNL/sundials/releases/download/v${PACKAGEVERSION}/sundials-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean
