################################################################
####
#### Makefile for Dealii installation
####
################################################################

PACKAGE = DEALII
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 9.1.1
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = mpi
MODULES = petsc

CONFIGURATION = default
include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D DEAL_II_WITH_MPI=ON \
    -D DEAL_II_CXX_FLAGS="-march=native" \
    -D DEAL_II_WITH_CXX17=OFF \
    -D DEAL_II_COMPONENT_EXAMPLES=OFF 
PKGCONFIGSET = lib/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

GITURL = https://github.com/dealii/dealii/archive/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.git

