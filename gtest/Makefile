################################################################
####
#### Makefile for googletest
####
################################################################

PACKAGE = googletest
URL = https://github.com/google/googletest
PACKAGEVERSION = 1.14.0
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build

TGZURL = https://github.com/google/googletest/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
