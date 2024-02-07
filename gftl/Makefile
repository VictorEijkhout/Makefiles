################################################################
####
#### Makefile for Mapl installation
####
################################################################

PACKAGE = gftl
URL = https://github.com/Goddard-Fortran-Ecosystem/gFTL
PACKAGEVERSION = 1.11.0
PACKAGEVERSION_MAJORMINOR = 1.11

MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
include ${MAKEINCLUDES}/Make.cmake
CMAKEPREFIXPATHSET = 1
INSTALLEXTRAPATH = GFTL-${PACKAGEVERSION_MAJORMINOR}
NOLIB = 1
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/Goddard-Fortran-Ecosystem/gFTL/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

