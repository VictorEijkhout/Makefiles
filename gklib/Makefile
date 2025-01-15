################################################################
####
#### Makefile for GKlib installation
####
#### https://github.com/scivision/METIS/issues/10
#### https://github.com/scivision/METIS/issues/8
####
################################################################

PACKAGE = gklib
URL = https://github.com/scivision/GKlib.git
PACKAGEVERSION = git
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
CMAKEPREFIXPATHSET = 1
CMAKEFLAGS =
# It says this in the github readme but does not work: --workflow --preset default

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

GITREPO = https://github.com/scivision/GKlib.git
GITDATE = today
include ${MAKEINCLUDES}/Make.git

