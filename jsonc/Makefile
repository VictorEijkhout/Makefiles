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

PKGCONFIG = lib64/pkgconfig
ifeq "${TACC_FAMILY_COMPILER}" "intel"
  CMAKECOMPILERFLAGS = -Wno-implicit-function-declaration
endif
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

GITREPO = https://github.com/json-c/json-c.git
include ${MAKEINCLUDES}/Make.git

