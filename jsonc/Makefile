################################################################
####
#### Makefile for json-c installation
####
#### https://github.com/json-c/json-c
####
################################################################

PACKAGE = jsonc
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

PKGCONFIG = lib64/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

GITREPO = https://github.com/json-c/json-c.git
include ${MAKEINCLUDES}/Make.git

