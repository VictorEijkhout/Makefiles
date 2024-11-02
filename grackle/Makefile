################################################################
####
#### Makefile for Grackle installation
####
################################################################

PACKAGE = GRACKLE
URL = https://github.com/grackle-project/grackle
PACKAGEVERSION = git
# 3.3.0 does not have cmakelists. Hm.
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = 

PKGCONFIG = 
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

GITREPO = https://github.com/grackle-project/grackle.git
include ${MAKEINCLUDES}/Make.git

TGZURL = https://github.com/grackle-project/grackle/archive/refs/tags/grackle-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

