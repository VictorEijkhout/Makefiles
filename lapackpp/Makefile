################################################################
####
#### Makefile for Lapack++
####
################################################################

PACKAGE = LAPACKPP
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = git
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = seq
MODULES = blaspp

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake

include ${MAKEINCLUDES}/Make.cbuild

GITREPO= git@github.com:icl-utk-edu/lapackpp.git
include ${MAKEINCLUDES}/Make.git
