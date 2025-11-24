################################################################
####
#### Makefile for ROMS installation
####
#### does not compile with intel 24
#### support only for registered users
####
################################################################

PACKAGE = ROMS
ABOUT = Regional Ocean Modeling System
URL = https://github.com/myroms/roms
PACKAGEVERSION = 4.2
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

DOWNLOADURL = https://github.com/myroms/roms/archive/refs/tags/roms-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/myroms/roms.git
include ${MAKEINCLUDES}/Make.git
