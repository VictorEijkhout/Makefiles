################################################################
####
#### Makefile for Wfmash installation
####
#### THIS DOES NOT BUILD
####
################################################################

PACKAGE = wfmash
SKIPPACKAGE = 1
ABOUT = Wfmash
PACKAGEVERSION = 0.24.1
URL = https://github.com/waveygang/wfmash
MODE = seq
MODULE = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/waveygang/wfmash/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/waveygang/wfmash.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git
