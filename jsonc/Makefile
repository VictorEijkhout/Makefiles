################################################################
####
#### Makefile for json-c installation
####
################################################################

PACKAGE = jsonc
URL = https://github.com/json-c/json-c
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
ifeq "${TACC_FAMILY_COMPILER}" "intel"
  CMAKECOMPILERFLAGS = -Wno-implicit-function-declaration
endif
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build public

GITREPO = https://github.com/json-c/json-c.git
include ${MAKEINCLUDES}/Make.git

