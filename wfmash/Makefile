################################################################
####
#### Makefile for Wfmash installation
####
################################################################

PACKAGE = wfmash
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
include ${MAKEINCLUDES}/Make.vars
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

TGZURL = https://github.com/waveygang/wfmash/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

