################################################################
####
#### Makefile for GSLLITE installation
####
################################################################

PACKAGE = GSLLITE
PACKAGEVERSION = 0.42.0
URL = https://github.com/gsl-lite/gsl-lite
MODE = seq

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

TGZURL = https://github.com/gsl-lite/gsl-lite/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
