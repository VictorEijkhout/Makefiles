################################################################
####
#### Makefile for DEMANGLE installation
####
################################################################

PACKAGE = DEMANGLE
ABOUT = demangler
PACKAGEVERSION = 0.0.0
URL = https://github.com/llvm/llvm-project/tree/main/libcxxabi/src/demangle
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
