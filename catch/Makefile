################################################################
####
#### Makefile for Catch2 installation
####
################################################################

PACKAGE = CATCH2
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 3.1.1
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODULES = 

include ${MAKEINCLUDES}/Make.info

## PKGCONFIGSET = lib64/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/catchorg/Catch2/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/catchorg/Catch2.git
include ${MAKEINCLUDES}/Make.git
