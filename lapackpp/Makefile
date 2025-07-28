################################################################
####
#### Makefile for Lapack++
####
################################################################

PACKAGE = LAPACKPP
PACKAGEVERSION = git
MODE = seq
MODULES = blaspp

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

GITREPO= git@github.com:icl-utk-edu/lapackpp.git
include ${MAKEINCLUDES}/Make.git
