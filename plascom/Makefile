################################################################
####
#### Makefile for Plascom2 installation
####
#### only needed for the Plascom2 CSA project
####
################################################################

PACKAGE = plascom2
URL = https://github.com/chess-uiuc/plascom2
PACKAGEVERSION = git

MODE = mpi
MODULES = phdf5
# adios2 hypre metis mfem mutationpp precice petsc

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
CMAKEFLAGS = 

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

GITREPO = git@github.com:chess-uiuc/plascom2.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

