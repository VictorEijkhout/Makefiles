################################################################
####
#### Makefile for Peridigm installation
####
#### needs a modern enough C++: gcc12 or so.
####
################################################################

PACKAGE = PERIDIGM
URL = https://github.com/peridigm/peridigm
PACKAGEVERSION = git
MODE = mpi
MODULES = trilinos pnetcdf boost

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -D Trilinos_DIR:PATH=${TACC_TRILINOS_DIR} \
    -D CMAKE_CXX_FLAGS=\"-qopenmp -std=c++17\"

## -D CMAKE_CXX_FLAGS:STRING="-O2 -Wall -std=c++11 -pedantic -Wno-long-long -ftrapv -Wno-deprecated" 

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

GITREPO = https://github.com/peridigm/peridigm.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean

