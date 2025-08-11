################################################################
####
#### Makefile for pflogger installation
####
#### building block for MAPL
####
################################################################

PACKAGE = pflogger
ABOUT = A parallel Fortran logger
URL = https://github.com/Goddard-Fortran-Ecosystem/pFlogger
PACKAGEVERSION = 1.12.0
PACKAGEVERSION_MAJORMINOR = 1.12

MODE = mpi
MODULES = yafyaml

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
include ${MAKEINCLUDES}/Make.cmake
CMAKEPREFIXPATHSET = 1
INSTALLEXTRAPATH = PFLOGGER-${PACKAGEVERSION_MAJORMINOR}
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/Goddard-Fortran-Ecosystem/pFlogger/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

