################################################################
####
#### Makefile for Precice installation
####
#### build in jail fails:
#### https://github.com/precice/precice/issues/2303
####
################################################################

PACKAGE = precice
ABOUT = coupling library for partitioned multi-physics simulations
URL = https://github.com/precice/precice
PACKAGEVERSION = 3.2.0
# 3.3.0

MODE = mpi
MODULES = boost eigen petsc

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
CMAKEFLAGS = \
    -D BUILD_TESTING=OFF \
    -D PRECICE_FEATURE_MPI_COMMUNICATION=ON \
    -D PRECICE_FEATURE_PYTHON_ACTIONS=OFF \
    -D PRECICE_FEATURE_PETSC_MAPPING=ON

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

TGZURL = https://github.com/precice/precice/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = git@github.com:precice/precice.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

