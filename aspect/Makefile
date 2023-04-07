################################################################
####
#### Makefile for Aspect installation
####
#### https://aspect.geodynamics.org/
#### https://github.com/geodynamics/aspect
####
################################################################

PACKAGE = ASPECT
PACKAGEVERSION = 2.4.0
MODE = mpi
MODULES = dealii
# dealii newer than 9.3

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
	-D CMAKE_BUILD_TYPE=Release \
        -D DEAL_II_DIR=${LMOD_DEALII_DIR} \
	-Wno-dev
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

ASPECT_BIN = ${ASPECT_INSTALLATION}/bin

TGZURL = https://github.com/geodynamics/aspect/releases/download/v${PACKAGEVERSION}/aspect-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

