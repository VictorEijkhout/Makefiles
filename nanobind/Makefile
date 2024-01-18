################################################################
####
#### Makefile for NANOBIND installation
####
################################################################

PACKAGE = NANOBIND
URL = https://github.com/wjakob/nanobind \
    https://nanobind.readthedocs.io/en/latest/
PACKAGEVERSION = git
# 1.8.0
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
CMAKEFLAGS = 
NOLIB = 1
INCLUDELOC = share/nanobind/include
CMAKEPREFIXPATHSET = 1
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build public
info ::
	@echo "suggestion to install as submodule"
	@echo "    git submodule add https://github.com/wjakob/nanobind ext/nanobind"
	@echo "    git submodule update --init --recursive"

GITREPO = https://github.com/wjakob/nanobind.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git

TGZURL = https://github.com/wjakob/nanobind/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
