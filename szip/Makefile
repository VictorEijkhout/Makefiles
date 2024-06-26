################################################################
####
#### Makefile for SZIP installation
####
################################################################

PACKAGE = SZIP
PACKAGEVERSION = 2.1.1
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKFLAGS=\
  -D CMAKE_C_FLAGS=-fPIC \
  -D CMAKE_CXX_FLAGS=-fPIC
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
CMAKEPREFIXPATHSET = 1
CMAKEMODULEPATHSET = share/cmake
default_install : configure build

TGZURL = https://docs.hdfgroup.org/archive/support/ftp/lib-external/szip/${PACKAGEVERSION}/src/szip-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

# there is something wrong with this installation:
# CMake Error at /Users/eijkhout/Installation/szip/build-2.1.1-macbookair-gcc14/szip-config.cmake:61 (include):
#   include could not find requested file:

#     /Users/eijkhout/Installation/share/cmake/szip-targets.cmake
