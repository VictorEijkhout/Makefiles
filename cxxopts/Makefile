################################################################
####
#### Makefile for CXXOPTS installation
####
################################################################

PACKAGE = CXXOPTS
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 3.0.0
MODULES = 
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGSET = lib64/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/jarro2783/cxxopts/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
