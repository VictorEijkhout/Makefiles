################################################################
####
#### Makefile for Aspect installation
####
#### https://github.com/geodynamics/aspect
####
################################################################

PACKAGE = ASPECT
ABOUT = Advanced Solver for Planetary Evolution, Convection, and Tectonics
URL = https://aspect.geodynamics.org/ \
    https://github.com/geodynamics/aspect
PACKAGEVERSION = 3.0.0
MODE = mpi
MODULES = dealii
# dealii newer than 9.3

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
CMAKEFLAGS = \
    -D CMAKE_EXE_LINKER_FLAGS=-Wl,-demangle \
    -D DEAL_II_DIR=${TACC_DEALII_DIR}
CMAKEBUILDTYPE = DebugRelease
include ${MAKEINCLUDES}/Make.cbuild
CMAKEPREFIXPATHSET = 1
CPTOINSTALLDIR = cookbooks
DEPENDSONCURRENT = boost

# libWorldBuilder is in BIN
# https://github.com/geodynamics/aspect/issues/6312

info ::
	@echo "make default_install : plain configure build"
.PHONY: default_install
default_install : configure build

HASBIN = 1

DOWNLOADURL = https://github.com/geodynamics/aspect/releases/download/v${PACKAGEVERSION}/aspect-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
include ${MAKEINCLUDES}/Make.clean

