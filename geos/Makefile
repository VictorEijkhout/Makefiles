################################################################
####
#### Makefile for Geos installation
####
#### installation needs a shitload of dependencies
####
################################################################

PACKAGE = GEOS
URL = https://github.com/GEOS-DEV/GEOS \
    https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/buildGuide/Index.html
PACKAGEVERSION = 1.2.0
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = 
CMAKESOURCE = src/main
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

DOWNLOADURL = https://github.com/GEOS-DEV/GEOS/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/GEOS-DEV/GEOS.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean

