################################################################
####
#### Makefile for Mapl installation
####
################################################################

PACKAGE = mapl
URL = https://github.com/GEOS-ESM/MAPL
PACKAGEVERSION = 2.43.2

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

##
## cmake version
##
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/GEOS-ESM/MAPL/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

