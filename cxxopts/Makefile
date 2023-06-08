################################################################
####
#### Makefile for CXXOPTS installation
####
################################################################

PACKAGE = CXXOPTS
PACKAGEVERSION = 3.1.1
MODE = seq

include ${MAKEINCLUDES}/Make.info

ifeq "${TACC_FAMILY_COMPILER}" "intel"
  LIBDIR = lib64
  PKGCONFIG = lib64/pkgconfig
else
  LIBDIR = lib
  PKGCONFIG = lib/pkgconfig
endif
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/jarro2783/cxxopts/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
