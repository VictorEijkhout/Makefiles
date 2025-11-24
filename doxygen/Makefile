################################################################
####
#### Makefile for doxygen installation
####
################################################################

PACKAGE = doxygen
ABOUT = Widely-used documentation generator tool in software development.
PACKAGEVERSION = 1.12.0
# 1.14.0 is out
PACKAGEVERSIONDOWNLOAD = 1_12_0
URL = https://github.com/doxygen/doxygen
MODE = core

include ${MAKEINCLUDES}/Make.info

CMAKEPREFIXPATHSET = 1
PKGCONFIGLIB = pkgconfig
HASBIN = 1
CMAKEFLAGS = \
    -D SPDLOG_INSTALL=ON
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build 

DOWNLOADURL = https://github.com/doxygen/doxygen/archive/refs/tags/Release_${PACKAGEVERSIONDOWNLOAD}.tar.gz
include ${MAKEINCLUDES}/Make.download
