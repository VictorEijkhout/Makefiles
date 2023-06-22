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
include ${MAKEINCLUDES}/Make.cmake

include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://download.osgeo.org/proj/proj-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
