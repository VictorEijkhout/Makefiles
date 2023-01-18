################################################################
####
#### Makefile for Dealii installation
####
################################################################

PACKAGE = DEALII
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 9.4.1
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = mpi
MODULES = p4est petsc trilinos

CONFIGURATION = default
include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D DEAL_II_WITH_MPI=ON \
    -D DEAL_II_CXX_FLAGS="-march=native" \
    -D DEAL_II_WITH_CXX17=OFF \
    -D DEAL_II_COMPONENT_EXAMPLES=OFF \
    -D DEAL_II_WITH_P4EST=ON \
    -D DEAL_II_WITH_PETSC=ON -D P4EST_DIR=${LMOD_P4EST_DIR} \
    -D DEAL_II_WITH_TRILINOS=ON -D TRILINOS_DIR=${LMOD_TRILINOS_DIR}
# -D DEAL_II_WITH_SLEPC=ON
PKGCONFIGSET = lib/pkgconfig
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: complex
info ::
	@echo "================"
	@echo "make complex: configure build with complex. Make sure to load complex petsc!"
complex :
	@make --no-print-directory configure build PACKAGEVERSION=${PACKAGEVERSION} \
	    INSTALLEXT="complex" CMAKEFLAGS="${CMAKEFLAGS} -DDEAL_II_WITH_PETSC_COMPLEX=ON"

GITREPO = https://github.com/dealii/dealii/archive/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.git

TGZURL = https://github.com/dealii/dealii/releases/download/v${PACKAGEVERSION}/dealii-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
