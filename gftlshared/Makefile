################################################################
####
#### Makefile for gFTL-shared installation
####
#### building block for MAPL
####
################################################################

PACKAGE = gftl-shared
URL = https://github.com/Goddard-Fortran-Ecosystem/gFTL-shared
PACKAGEVERSION = 1.7.0
PACKAGEVERSION_MAJORMINOR = 1.7

MODE = seq
MODULES = gftl

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
include ${MAKEINCLUDES}/Make.cmake
MODULENAME = gftlshared
CMAKEPREFIXPATHSET = 1
INSTALLEXTRAPATH = GFTL_SHARED-${PACKAGEVERSION_MAJORMINOR}
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
default_install : configure build

DOWNLOADURL = https://github.com/Goddard-Fortran-Ecosystem/gFTL-shared/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

