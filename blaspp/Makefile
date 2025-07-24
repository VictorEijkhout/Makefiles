################################################################
####
#### Makefile for Blas++
####
################################################################

PACKAGE = BLASPP
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
CMAKEFLAGS = -Dblas=mkl
LIBDIR = lib64
# .cmake files
CMAKEPREFIXPATHSET = lib64/blaspp

include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

GITREPO= git@github.com:icl-utk-edu/blaspp.git
include ${MAKEINCLUDES}/Make.git
