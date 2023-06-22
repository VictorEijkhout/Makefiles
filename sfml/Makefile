################################################################
####
#### Makefile for SFML installation
####
################################################################

PACKAGE = SFML
URL = https://github.com/SFML/SFML
PACKAGEVERSION = 2.6.0
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
## CMAKEFLAGS = -D CMAKE_COMPILE_WARNING_AS_ERROR=OFF
## CMAKECOMPILERFLAGS = -Wno-error
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/SFML/SFML/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
