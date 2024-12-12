################################################################
####
#### Makefile for Microsoft Guideline Suport Library
####
################################################################

PACKAGE = msgsl
PACKAGEVERSION = 4.1.0
URL = https://github.com/microsoft/GSL
MODE = seq

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = -D GSL_TEST=OFF
include ${MAKEINCLUDES}/Make.cmake
NOLIB = 1
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/microsoft/GSL/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
