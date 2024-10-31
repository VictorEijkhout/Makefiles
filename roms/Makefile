################################################################
####
#### Makefile for ROMS installation
####
################################################################

PACKAGE = ROMS
URL = https://roms.org/index.html
PACKAGEVERSION = 4.1
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
HASBIN = 1
CMAKEFLAGS = 
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build public

TGZURL = https://github.com/myroms/roms/archive/refs/tags/roms-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/myroms/roms.git
include ${MAKEINCLUDES}/Make.git
