################################################################
####
#### Makefile for Mutationpp installation
####
#### only needed for the CHyPS CSA project
####
################################################################

PACKAGE = mutationpp
ABOUT = The MUlticomponent Thermodynamic And Transport library for IONized gases in C++
URL = https://github.com/mutationpp/Mutationpp
PACKAGEVERSION = git

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

##
## cmake 
##
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

GITREPO = git@github.com:mutationpp/Mutationpp.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

