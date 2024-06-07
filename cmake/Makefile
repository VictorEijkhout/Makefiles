################################################################
####
#### Makefile for CMAKE installation
####
################################################################

PACKAGE = CMAKE
PACKAGEVERSION = 3.29.5
URL = https://cmake.org/
MODE = core

include ${MAKEINCLUDES}/Make.info

CMAKEPREFIXPATHSET = 1
NOLIB = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build 

TGZURL = https://github.com/Kitware/CMake/releases/download/v${PACKAGEVERSION}/cmake-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
