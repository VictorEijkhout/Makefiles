################################################################
####
#### Makefile for Libsc installation
####
#### This is needed for p4est, but should be included
####
################################################################

PACKAGE = libsc
ABOUT = The "sc" auxiliary library
URL = https://github.com/cburstedde/libsc
PACKAGEVERSION = 2.8.7
SKIPPACKAGE = 1

MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

##
## Cmake
##
CMAKEFLAGS =
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

##
## download
##
TGZURL = https://github.com/cburstedde/libsc/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
