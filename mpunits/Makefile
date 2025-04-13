################################################################
####
#### Makefile for MPUNITS installation
####
################################################################

PACKAGE = MPUNITS
PACKAGEVERSION = 2.4.0
URL = https://github.com/mpusz/mp-units
MODE = seq
MODULES = catch2 fmtlib gsllite

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
NOLIB = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build public

TGZURL = https://github.com/mpusz/mp-units/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

