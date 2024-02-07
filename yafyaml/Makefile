################################################################
####
#### Makefile for yafyaml installation
####
#### building block for MAPL
####
################################################################

PACKAGE = yafyaml
URL = https://github.com/Goddard-Fortran-Ecosystem/yaFyaml
PACKAGEVERSION = 1.2.0
PACKAGEVERSION_MAJORMINOR = 1.2

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
include ${MAKEINCLUDES}/Make.cmake
CMAKEPREFIXPATHSET = 1
INSTALLEXTRAPATH = YAFYAML-${PACKAGEVERSION_MAJORMINOR}
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/Goddard-Fortran-Ecosystem/yaFyaml/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

