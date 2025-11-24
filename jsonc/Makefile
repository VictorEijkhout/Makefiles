################################################################
####
#### Makefile for json-c installation
####
################################################################

PACKAGE = jsonc
ABOUT = A JSON implementation in C
URL = https://github.com/json-c/json-c
PACKAGEVERSION = 0.18
PACKAGEVERSION_DOWNLOAD = 0.18-20240915
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1

CMAKEFLAGS = \
    -D DISABLE_WERROR:BOOL=ON

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
default_install : configure build public

DOWNLOADURL = https://github.com/json-c/json-c/archive/refs/tags/json-c-${PACKAGEVERSION_DOWNLOAD}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/json-c/json-c.git
include ${MAKEINCLUDES}/Make.git

