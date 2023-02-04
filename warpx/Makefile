################################################################
####
#### Makefile for WarpX
####
################################################################

PACKAGE = warpx
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = git
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = mpi
MODULES = blaspp lapackpp

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake

include ${MAKEINCLUDES}/Make.cbuild

GITREPO= https://github.com/ECP-WarpX/WarpX.git
include ${MAKEINCLUDES}/Make.git
