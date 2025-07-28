################################################################
####
#### Makefile for GMP installation
####
################################################################

PACKAGE = GMP
URL = https://gmplib.org/
PACKAGEVERSION = 6.3.0
MODE = core
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
CMAKEPREFIXPATHSET = 1
PKGCONFIGLIB = pkgconfig

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build 

TXZURL = https://gmplib.org/download/gmp/gmp-${PACKAGEVERSION}.tar.xz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
