################################################################
####
#### Makefile for CXXOPTS installation
####
################################################################

PACKAGE = CXXOPTS
PACKAGEVERSION = 3.2.0
URL = https://github.com/jarro2783/cxxopts
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

TGZURL = https://github.com/jarro2783/cxxopts/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/jarro2783/cxxopts.git
