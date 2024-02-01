################################################################
####
#### Makefile for SZLIB installation
####
################################################################

PACKAGE = SZLIB
PACKAGEVERSION = 2.1.1
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKFLAGS=\
  -D CMAKE_C_FLAGS=-fPIC \
  -D CMAKE_CXX_FLAGS=-fPIC
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build

TGZURL = https://docs.hdfgroup.org/archive/support/ftp/lib-external/szip/${PACKAGEVERSION}/src/szip-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean
