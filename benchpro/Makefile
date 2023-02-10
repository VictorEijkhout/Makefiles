################################################################
####
#### Makefile for benchpro
####
################################################################

PACKAGE = benchpro
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = git
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = seq

include ${MAKEINCLUDES}/Make.info

configure build :
	@echo "no such rule"
include ${MAKEINCLUDES}/Make.vars

GITREPO = git clone https://github.com/TACC/benchpro.git
include ${MAKEINCLUDES}/Make.git
