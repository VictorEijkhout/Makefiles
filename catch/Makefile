################################################################
####
#### Makefile for Catch2 installation
####
################################################################

PACKAGE = CATCH2
PACKAGEVERSION = 3.1.1
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

## PKGCONFIGSET = lib64/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/catchorg/Catch2/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/catchorg/Catch2.git
include ${MAKEINCLUDES}/Make.git
