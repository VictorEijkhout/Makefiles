################################################################
####
#### Makefile for Mutationpp installation
####
#### only needed for the CHyPS CSA project
####
################################################################

PACKAGE = mutationpp
URL = https://github.com/mutationpp/Mutationpp
PACKAGEVERSION = git

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
CMAKEFLAGS = -D CMAKE_BUILD_TYPE=Release
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

GITREPO = git@github.com:mutationpp/Mutationpp.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

