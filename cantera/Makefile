################################################################
####
#### Makefile for CANTERA installation
####
################################################################

PACKAGE = CANTERA
ABOUT = 
PACKAGEVERSION = 3.1.0
URL = https://github.com/Cantera/cantera
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

####
#### Download
####
TGZURL = https://github.com/Cantera/cantera/releases/download/v${PACKAGEVERSION}/cantera-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/Cantera/cantera.git
include ${MAKEINCLUDES}/Make.git
