################################################################
####
#### Makefile for fmtlib
####
################################################################

PACKAGE = CXXOPTS
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 6.1.2
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODULES = 
MODULEROOT =${HOME}/modulefiles

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/fmtlib/fmt/archive/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
