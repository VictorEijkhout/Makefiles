################################################################
####
#### Makefile for ZLIB installation
####
################################################################

PACKAGE = ZLIB
URL = https://github.com/madler/zlib
PACKAGEVERSION = 1.3.1
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKFLAGS=\
  -D CMAKE_C_FLAGS=-fPIC \
  -D CMAKE_CXX_FLAGS=-fPIC
include ${MAKEINCLUDES}/Make.cmake
PKGCONFIG = share/pkgconfig
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/madler/zlib/releases/download/v${PACKAGEVERSION}/zlib-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean
