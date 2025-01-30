################################################################
####
#### Makefile for Precice installation
####
#### only needed for the CHyPS CSA project
####
################################################################

PACKAGE = precice
URL = https://github.com/precice/precice
PACKAGEVERSION = 3.1.2
# 2.5.0

MODE = mpi
MODULES = boost eigen petsc

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
CMAKEFLAGS = \
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

