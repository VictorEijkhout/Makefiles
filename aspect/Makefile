################################################################
####
#### Makefile for Aspect installation
####
################################################################

PACKAGE = ASPECT
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 2.2.0
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = mpi
MODULES = dealii

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

ASPECT_BIN = ${ASPECT_INSTALLATION}/bin

CMAKEFLAGS = \
        -D DEAL_II_DIR=${TACC_DEALII_DIR} \

include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

