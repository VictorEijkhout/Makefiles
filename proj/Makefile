################################################################
####
#### Makefile for PROJ installation
#### https://proj.org/index.html
####
################################################################

PACKAGE = PROJ
PACKAGEVERSION = 9.2.1
MODE = seq
MODULES = sqlite

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = -D BUILD_TESTING=OFF
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cmake
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build public

TGZURL = https://download.osgeo.org/proj/proj-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
