################################################################
####
#### Makefile for fArgParse installation
####
#### building block for MAPL
####
################################################################

PACKAGE = fargparse
URL = https://github.com/Goddard-Fortran-Ecosystem/fArgParse
PACKAGEVERSION = 1.6.0
PACKAGEVERSION_MAJORMINOR = 1.6

MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
include ${MAKEINCLUDES}/Make.cmake
CMAKEPREFIXPATHSET = 1
INSTALLEXTRAPATH = FARGPARSE-${PACKAGEVERSION_MAJORMINOR}
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
default_install : configure build

DOWNLOADURL = https://github.com/Goddard-Fortran-Ecosystem/fArgParse/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

