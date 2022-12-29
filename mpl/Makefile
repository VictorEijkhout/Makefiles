################################################################
####
#### Makefile for mpl
####
################################################################

PACKAGE = MPL
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 0.2.0
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

GITREPO = git@github.com:rabauke/mpl.git
include ${MAKEINCLUDES}/Make.git

TGZURL = https://github.com/rabauke/mpl/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
