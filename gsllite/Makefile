################################################################
####
#### Makefile for GSLLITE installation
####
################################################################

PACKAGE = GSLLITE
PACKAGEVERSION = 0.43.0
# 1.0.1 is out
URL = https://github.com/gsl-lite/gsl-lite
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
NOLIB = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build public

TGZURL = https://github.com/gsl-lite/gsl-lite/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
