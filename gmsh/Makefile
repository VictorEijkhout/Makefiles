################################################################
####
#### Makefile for GMSH installation
####
################################################################

PACKAGE = GMSH
ABOUT = Mesh Generator
PACKAGEVERSION = 4.15.0
URL = https://gmsh.info/
MODE = omp

include ${MAKEINCLUDES}/Make.info

NOINC = 1
NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build public

DOWNLOADURL = https://gmsh.info/src/gmsh-${PACKAGEVERSION}-source.tgz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
