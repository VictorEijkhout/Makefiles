################################################################
####
#### Makefile for Aspect installation
####
#### https://github.com/geodynamics/aspect
####
################################################################

PACKAGE = ASPECT
URL = https://aspect.geodynamics.org/ \
    https://github.com/geodynamics/aspect
PACKAGEVERSION = 3.0.0
#2.5.0
MODE = mpi
MODULES = dealii
# dealii newer than 9.3

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
        -D DEAL_II_DIR=${TACC_DEALII_DIR}
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
info ::
	@echo "make default_install : plain configure build"
default_install : configure build

HASBIN = 1

TGZURL = https://github.com/geodynamics/aspect/releases/download/v${PACKAGEVERSION}/aspect-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

