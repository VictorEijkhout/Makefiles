################################################################
####
#### Makefile for Eigen installation
####
################################################################

PACKAGE = ninja
PACKAGEVERSION = 1.12.1
URL = https://github.com/ninja-build/ninja
MODULES = 
MODE = core

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
HASBIN = 1
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cbuild

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build public

TGZURL = https://github.com/ninja-build/ninja/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
