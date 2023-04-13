################################################################
####
#### Makefile for Peridigm installation
####
#### https://peridigm.geodynamics.org/
#### https://github.com/geodynamics/peridigm
####
################################################################

PACKAGE = PERIDIGM
PACKAGEVERSION = git
MODE = mpi
MODULES = trilinos boost

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
-D CMAKE_BUILD_TYPE:STRING=Release \
-D Trilinos_DIR:PATH=${TACC_TRILINOS_DIR} \
-D BOOST_ROOT=${TACC_BOOST_DIR} \
-D CMAKE_CXX_FLAGS=\"-Wno-dev -qopenmp -std=c++17\"

## -D CMAKE_CXX_FLAGS:STRING="-O2 -Wall -std=c++11 -pedantic -Wno-long-long -ftrapv -Wno-deprecated" 

include ${MAKEINCLUDES}/Make.cmake
HASBIN=1
include ${MAKEINCLUDES}/Make.cbuild

GITREPO = https://github.com/peridigm/peridigm.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean

