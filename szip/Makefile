################################################################
####
#### Makefile for ZLIB installation
####
################################################################

PACKAGE = ZLIB
PACKAGEVERSION = 1.2.13
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKFLAGS=\
  -D CMAKE_C_FLAGS=-fPIC \
  -D CMAKE_CXX_FLAGS=-fPIC
include ${MAKEINCLUDES}/Make.cmake
PKGCONFIG = share/pkgconfig
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/madler/zlib/releases/download/v${PACKAGEVERSION}/zlib-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean
