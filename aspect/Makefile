################################################################
####
#### Makefile for Aspect installation
####
#### https://aspect.geodynamics.org/
####
################################################################

PACKAGE = ASPECT
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 2.4.0
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = mpi
MODULES = dealii

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
	-D CMAKE_BUILD_TYPE=Release \
        -D DEAL_II_DIR=${TACC_DEALII_DIR} \
	-Wno-dev
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

ASPECT_BIN = ${ASPECT_INSTALLATION}/bin

TGZURL = https://github.com/geodynamics/aspect/releases/download/v${PACKAGEVERSION}/aspect-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

