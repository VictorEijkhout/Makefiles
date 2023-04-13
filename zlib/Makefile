################################################################
####
#### Makefile for ZLIB installation
####
################################################################

PACKAGE = ZLIB
PACKAGEVERSION = 1.2.13
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGSET = share/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/madler/zlib/releases/download/v${PACKAGEVERSION}/zlib-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean
